#!/usr/bin/env bash
set -euo pipefail

############################################################
# Bulk directory-only ACL setter for Globus (GCS v5) via CLI
#
# - Enumerates directories (only) beneath a base path
# - Applies ACLs for multiple principals (users or groups)
# - Excludes specified subfolders; supports regex include/exclude
# - Depth control (e.g., -d 1 for first-level only)
# - Dry-run and verbose modes
# - Portable: no 'mapfile'; no associative arrays (Bash 3+ OK)
#
# Pre-flight uses: 'globus endpoint permission list' (ACL subsystem)
#
# Requires: Globus CLI ('globus'), jq, bash
# Login first: 'globus login'
############################################################

# Defaults
ENDPOINT_ID=""
BASE_PATH=""
PERMISSION="rw"                  # allowed: r | rw
PRINCIPAL_TYPE="identity"        # identity | group
USERS_FILE=""
USERS_CSV=""
EXCLUDE_DIRS=()                  # comma-separated names under BASE_PATH to exclude
REGEX_INCLUDE=".*"               # include all by default
REGEX_EXCLUDE=""                 # none by default
MAX_DEPTH=-1                     # -1: unlimited; 1: first level; 2, ...
DRY_RUN=false
VERBOSE=false

print_usage() {
  cat <<'EOF'
Usage:
  globus_acl_bulk.sh -e COLLECTION_ID -b BASE_PATH (-U users.txt | -u "a,b,c") [options]

Required:
  -e COLLECTION_ID      Globus GCS v5 Guest/Mapped Collection UUID (e.g., 0c239687-345e-420e-9d1d-f0002604658f)
  -b BASE_PATH          Base folder (e.g., /CM)
  -U USERS_FILE         File with one principal per line (email or identity UUID, or group UUID)
  -u USERS_CSV          Comma-separated principals (emails/UUIDs). May be combined with -U.

Options:
  -p PERMISSION         Permissions: r | rw  (default: rw)
  -t TYPE               Principal type: identity | group  (default: identity)
  -x EXCLUDE_LIST       Comma-separated subfolder names under BASE_PATH to exclude (e.g., "x,y")
  -I REGEX_INCLUDE      Regex applied to FULL absolute paths to include (default: ".*")
  -X REGEX_EXCLUDE      Regex applied to FULL absolute paths to exclude
  -d MAX_DEPTH          Max depth from BASE_PATH: -1 (unlimited), 1 (first level), 2, ... (default: -1)
  -n                    Dry run (print planned ACLs; no changes)
  -v                    Verbose logging
  -h                    Help

Examples:
  # RW for users A,B on ALL dirs under /CM except 'pwet' and 'zbla' (first level only)
  ./globus_acl.sh -e <COLL_ID> -b /CM -u "a@ex.com,b@ex.com" -p rw -x "pwet,zbla" -d 1 -v

  # Read-only for a group on /CM (all depths), include only cohorts 'pwet' and 'zbla', exclude '/CM/truite'
  ./globus_acl.sh -e <COLL_ID> -b /CM -u "GROUP_UUID" -t group -p r \
     -I "/CM/(pwet|zbla)/" -X "/CM/truite/" -v
EOF
}

log()  { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
vlog() { $VERBOSE && log "$@" || true; }
die()  { echo "ERROR: $*" >&2; exit 1; }

# ---- Parse args
while getopts ":e:b:U:u:p:t:x:I:X:d:nvh" opt; do
  case "$opt" in
    e) ENDPOINT_ID="$OPTARG" ;;
    b) BASE_PATH="$OPTARG" ;;
    U) USERS_FILE="$OPTARG" ;;
    u) USERS_CSV="$OPTARG" ;;
    p) PERMISSION="$OPTARG" ;;
    t) PRINCIPAL_TYPE="$OPTARG" ;;
    x) IFS=',' read -r -a EXCLUDE_DIRS <<< "$OPTARG" ;;
    I) REGEX_INCLUDE="$OPTARG" ;;
    X) REGEX_EXCLUDE="$OPTARG" ;;
    d) MAX_DEPTH="$OPTARG" ;;
    n) DRY_RUN=true ;;
    v) VERBOSE=true ;;
    h) print_usage; exit 0 ;;
    \?) die "Invalid option: -$OPTARG" ;;
    :)  die "Option -$OPTARG requires an argument." ;;
  esac
done

# ---- Validate args
[[ -z "$ENDPOINT_ID" || -z "$BASE_PATH" ]] && { print_usage; die "-e and -b are required."; }
if [[ -z "$USERS_FILE" && -z "$USERS_CSV" ]]; then
  print_usage; die "Provide principals via -U users.txt or -u \"a,b,c\".";
fi
[[ -n "$USERS_FILE" && ! -f "$USERS_FILE" ]] && die "Users file not found: $USERS_FILE"
[[ "$PRINCIPAL_TYPE" != "identity" && "$PRINCIPAL_TYPE" != "group" ]] && die "-t must be identity|group"
case "$PERMISSION" in r|rw) ;; *) die "Globus ACLs accept only -p r|rw";; esac
[[ "$MAX_DEPTH" =~ ^-?[0-9]+$ ]] || die "-d MAX_DEPTH must be an integer (e.g., 1, 2, or -1)."

# Normalize base path (no trailing slash)
BASE_PATH="${BASE_PATH%/}"

log "Collection: $ENDPOINT_ID"
log "Base path: $BASE_PATH"
log "Permission: $PERMISSION"
log "Principal type: $PRINCIPAL_TYPE"
$VERBOSE && log "Exclude subfolders under base: ${EXCLUDE_DIRS[*]:-<none>}"
$VERBOSE && log "Regex include: $REGEX_INCLUDE"
$VERBOSE && log "Regex exclude: ${REGEX_EXCLUDE:-<none>}"
$VERBOSE && log "Max depth: $MAX_DEPTH"
$DRY_RUN && log "Dry-run mode enabled (no changes will be made)."

# ---- Pre-flight: Globus CLI + jq + login
if ! command -v globus >/dev/null 2>&1; then
  die "Globus CLI not found. Please install Globus CLI and run 'globus login'. (On HPCs, try: module load Globus-CLI)"
fi
if ! command -v jq >/dev/null 2>&1; then
  die "'jq' not found. Please install jq (or load it via your environment modules)."
fi
if ! globus whoami >/dev/null 2>&1; then
  die "You are not logged in to Globus. Run 'globus login' before using this script."
fi

# ---- Pre-flight: ACL subsystem reachable (authoritative for your environment)
if ! globus endpoint permission list "$ENDPOINT_ID" --format json >/dev/null 2>&1; then
  die "ACL subsystem not accessible for this collection. Check session consent or your privileges."
fi
vlog "ACL subsystem reachable via 'globus endpoint permission list'."


# ---- Cache all existing ACLs once (JSON)
ALL_ACLS_JSON="$(globus endpoint permission list "$ENDPOINT_ID" --format json 2>/dev/null || echo '{}')"

# ---- Helper: normalize directory path (ensure trailing '/')
norm_dir_path() {
  local p="$1"
  echo "${p%/}/"
}

# ---- Helper: test if an exact ACL already exists
# Args: endpoint_id path principal_type principal_uuid permission
acl_exists_exact() {
  local ep="$1" path="$2" ptype="$3" puuid="$4" perm="$5"
  local path_norm
  path_norm="$(norm_dir_path "$path")"
  echo "$ALL_ACLS_JSON" \
    | jq -e --arg path "$path_norm" --arg ptype "$ptype" --arg puuid "$puuid" --arg perm "$perm" '
        .DATA[]?
        | select(.path == $path
                 and .principal_type == $ptype
                 and .principal == $puuid
                 and .permissions == $perm)
      ' >/dev/null
}

# ---- Helper: test if a superset ACL exists (rw satisfies r)
# Returns success if requested 'r' and existing rule is 'rw' on same path and principal
acl_exists_superset() {
  local ep="$1" path="$2" ptype="$3" puuid="$4" perm="$5"
  [[ "$perm" != "r" ]] && return 1  # only applicable when requesting 'r'
  local path_norm
  path_norm="$(norm_dir_path "$path")"
  echo "$ALL_ACLS_JSON" \
    | jq -e --arg path "$path_norm" --arg ptype "$ptype" --arg puuid "$puuid" '
        .DATA[]?
        | select(.path == $path
                 and .principal_type == $ptype
                 and .principal == $puuid
                 and .permissions == "rw")
      ' >/dev/null
}

# ---- Load principals
PRINCIPALS=()
if [[ -n "$USERS_FILE" ]]; then
  while IFS= read -r u || [[ -n "$u" ]]; do
    [[ -z "$u" || "$u" =~ ^# ]] && continue
    PRINCIPALS+=( "$u" )
  done < "$USERS_FILE"
fi
if [[ -n "$USERS_CSV" ]]; then
  IFS=',' read -r -a CSV_ARRAY <<< "$USERS_CSV"
  for u in "${CSV_ARRAY[@]}"; do
    u_trim="${u//[[:space:]]/}"
    [[ -z "$u_trim" ]] && continue
    PRINCIPALS+=( "$u_trim" )
  done
fi
[[ ${#PRINCIPALS[@]} -eq 0 ]] && die "No principals parsed from -U/-u input."
$VERBOSE && { log "Principals:"; for p in "${PRINCIPALS[@]}"; do echo "  - $p"; done; }

# ---- Validate principals and print names
PRINCIPAL_IDS=()
for USER in "${PRINCIPALS[@]}"; do
  if [[ "$PRINCIPAL_TYPE" == "identity" ]]; then
    ID_JSON="$(globus get-identities "$USER" --format json 2>/dev/null || true)"
    ID="$(echo "$ID_JSON" | jq -r '.identities[0].id // empty')"
    FULLNAME="$(echo "$ID_JSON" | jq -r '.identities[0].name // empty')"
    USERNAME="$(echo "$ID_JSON" | jq -r '.identities[0].username // empty')"
    [[ -z "$ID" ]] && die "Identity not found or invalid: $USER"
    [[ -z "$FULLNAME" ]] && FULLNAME="<no name available>"
    PRINCIPAL_IDS+=( "$ID" )
    log "Validated identity: $USER"
    log "  → UUID:      $ID"
    log "  → Username:  $USERNAME"
    log "  → Full Name: $FULLNAME"
  else
    # group
    if ! GROUP_JSON="$(globus group show "$USER" --format json 2>/dev/null)"; then
      die "Group not found or invalid: $USER"
    fi
    GROUP_NAME="$(echo "$GROUP_JSON" | jq -r '.name // empty')"
    [[ -z "$GROUP_NAME" ]] && GROUP_NAME="<unnamed group>"
    PRINCIPAL_IDS+=( "$USER" )  # group input is already a UUID
    log "Validated group: $USER"
    log "  → Group Name: $GROUP_NAME"
  fi
done

# ---- Helper: list immediate child directory names under a given absolute path (JSON + jq)
list_child_dir_names() {
  local parent="$1"
  globus ls "${ENDPOINT_ID}:${parent}/" --format json \
    | jq -r '.DATA[] | select(.type=="dir") | .name'
}

# ---- Build list of absolute excluded subpaths
EXCLUDE_ABS=()
for d in "${EXCLUDE_DIRS[@]}"; do
  [[ -z "$d" ]] && continue
  EXCLUDE_ABS+=( "${BASE_PATH%/}/$d" )
done

# ---- Traverse breadth-first; collect ALL_DIRS (portable queue)
QUEUE_PATHS=( "$BASE_PATH" )
QUEUE_DEPTHS=( 0 )
ALL_DIRS=()

queue_len() { echo "${#QUEUE_PATHS[@]}"; }

while (( $(queue_len) > 0 )); do
  # pop front
  cur="${QUEUE_PATHS[0]}"
  depth="${QUEUE_DEPTHS[0]}"
  QUEUE_PATHS=( "${QUEUE_PATHS[@]:1}" )
  QUEUE_DEPTHS=( "${QUEUE_DEPTHS[@]:1}" )

  # Stop descending if depth limit reached
  if (( MAX_DEPTH >= 0 && depth >= MAX_DEPTH )); then
    vlog "Depth $depth reached at $cur (limit $MAX_DEPTH) — not descending"
    continue
  fi

  # Pull child dir names via JSON + jq
  CHILD_NAMES_STR="$(list_child_dir_names "$cur" || true)"

  while IFS= read -r name; do
    [[ -z "$name" ]] && continue
    child="$cur/${name%/}"
    child="${child%/}/"  # normalize to trailing slash


    # Skip excluded subtree prefixes
    skip=false
    for ex in "${EXCLUDE_ABS[@]}"; do
      if [[ -n "$ex" && "$child" == "$ex"* ]]; then
        vlog "Excluding subtree: $child (matches $ex)"
        skip=true
        break
      fi
    done
    $skip && continue

    # Keep child and enqueue for further descent
    ALL_DIRS+=( "$child" )
    QUEUE_PATHS+=( "$child" )
    QUEUE_DEPTHS+=( $((depth + 1)) )
  done <<< "$CHILD_NAMES_STR"
done

# ---- Apply include/exclude regex filters over FULL paths
FILTERED_DIRS=()
for dir in "${ALL_DIRS[@]}"; do
  if [[ "$dir" =~ $REGEX_INCLUDE ]]; then
    if [[ -n "$REGEX_EXCLUDE" && "$dir" =~ $REGEX_EXCLUDE ]]; then
      continue
    fi
    FILTERED_DIRS+=( "$dir" )
  fi
done
[[ ${#FILTERED_DIRS[@]} -eq 0 ]] && die "After filtering, no directories remain. Adjust -x/-I/-X/-d."

log "Candidate directories to receive ACLs: ${#FILTERED_DIRS[@]}"
$VERBOSE && printf '%s\n' "${FILTERED_DIRS[@]}" | sed 's/^/  - /'

# ---- Confirm and apply
log "Starting ACL operations in 3 seconds. Ctrl+C to abort."
sleep 3

ERRORS=0
for i in "${!PRINCIPALS[@]}"; do
  USER="${PRINCIPALS[$i]}"
  USER_ID="${PRINCIPAL_IDS[$i]}"
  log "Principal: $USER"
  for dir in "${FILTERED_DIRS[@]}"; do
    DIR_FOR_ACL="$(norm_dir_path "$dir")"
    TARGET="${ENDPOINT_ID}:${DIR_FOR_ACL}"
    # 1) Exact duplicate check
    if acl_exists_exact "$ENDPOINT_ID" "$DIR_FOR_ACL" "$PRINCIPAL_TYPE" "$USER_ID" "$PERMISSION"; then
      log "  → SKIP (already present): $PERMISSION for $PRINCIPAL_TYPE '$USER' on $DIR_FOR_ACL"
      continue
    fi
    # 2) Superset check (rw satisfies r)
    if acl_exists_superset "$ENDPOINT_ID" "$DIR_FOR_ACL" "$PRINCIPAL_TYPE" "$USER_ID" "$PERMISSION"; then
      log "  → SKIP (superset exists: rw): requested '$PERMISSION' for $PRINCIPAL_TYPE '$USER' on $DIR_FOR_ACL"
      continue
    fi
    # 3) Create if not present
    if $DRY_RUN; then
      vlog "DRY-RUN: globus endpoint permission create '$TARGET' --permissions '$PERMISSION' --$PRINCIPAL_TYPE '$USER'"
      continue
    fi
    $VERBOSE && log "Setting $PERMISSION for $PRINCIPAL_TYPE '$USER' on $DIR_FOR_ACL"
    if ! globus endpoint permission create \
           "$TARGET" \
           --permissions "$PERMISSION" \
           --"$PRINCIPAL_TYPE" "$USER"; then
      echo "ERROR: Failed ACL on $DIR_FOR_ACL (principal: $USER)" >&2
      ((ERRORS++))
    fi
  done
done

if (( ERRORS > 0 )); then
  echo "Completed with $ERRORS error(s)." >&2
  exit 2
fi

log "All ACL operations completed successfully."
