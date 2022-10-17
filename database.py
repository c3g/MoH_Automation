"""Module providing database tables and operations support."""
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    Boolean,
    String,
    JSON,
    DateTime
    )
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

Base = declarative_base()
# session = sessionmaker()(bind=engine)

class BaseTable(Base):
    """
    Define fields common of all tables in database
    BaseTable:
        id integer [PK]
        deleted boolean
        extra_metadata json
    """
    __abstract__ = True

    id = Column(Integer, primary_key=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    def __init__(self, deleted=False, extra_metadata=None):
        self.deleted = deleted
        self.extra_metadata = extra_metadata

class Project(BaseTable):
    """
    Project:
        id integer [PK]
        name text
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "project"

    name = Column(String, nullable=False, unique=True)

    def __init__(self, name=None):
        super().__init__()
        self.name = name

    def __repr__(self):
        return f"Project({self.name!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Patient(BaseTable):
    """
    Patient:
        id integer [PK]
        project_id integer [ref: > project.id]
        name text
        alias blob
        cohort text
        institution text
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "patient"

    project_id = Column(Integer, ForeignKey("project.id"))
    name = Column(String, nullable=False, unique=True)
    alias = Column(String, nullable=True)
    cohort = Column(String, nullable=True)
    institution = Column(String, nullable=True)

    project = relationship("Project", backref="patient", lazy=False)

    def __init__(self, name=None, alias=None, cohort=None, institution=None):
        super().__init__()
        self.name = name
        self.alias = alias
        self.cohort = cohort
        self.institution = institution

    def __repr__(self):
        return f"Patient({self.project!r}, {self.name!r}, {self.alias!r}, {self.cohort!r}, {self.institution!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Experiment(BaseTable):
    """
    Experiment:
        id integer [PK]
        run_id integer [ref: > run.id]
        project_id integer [ref: > project.id]
        sequencing_technology text
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "experiment"

    run_id = Column(Integer, ForeignKey("run.id"))
    project_id = Column(Integer, ForeignKey("project.id"))
    sequencing_technology = Column(String, nullable=True)

    run = relationship("Run", backref="experiment", lazy=False)
    project = relationship("Project", backref="experiment", lazy=False)

    def __init__(self, sequencing_technology=None):
        super().__init__()
        self.sequencing_technology = sequencing_technology

    def __repr__(self):
        return f"Experiment({self.run!r}, {self.project!r}, {self.sequencing_technology!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Sample(BaseTable):
    """
    Sample:
        id integer [PK]
        patient_id integer [ref: > patient.id]
        experiment_id integer [ref: > experiment.id]
        name text
        tumour boolean
        alias blob
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "sample"

    patient_id = Column(Integer, ForeignKey("patient.id"))
    experiment_id = Column(Integer, ForeignKey("experiment.id"))
    name = Column(String, nullable=False)
    tumour = Column(Boolean, default=False)
    alias = Column(String, nullable=True)

    patient = relationship("Patient", backref="sample", lazy=False)
    experiment = relationship("Experiment", backref="sample", lazy=False)

    def __init__(self, name=None, tumour=False, alias=None):
        super().__init__()
        self.name = name
        self.tumour = tumour
        self.alias = alias

    def __repr__(self):
        return f"Sample({self.patient!r}, {self.experiment!r}, {self.name!r}, {self.tumour!r}, {self.alias!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Run(BaseTable):
    """
    Patient:
        id integer [PK]
        lab_id text
        name text
        date timestamp
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "run"

    lab_id = Column(String, nullable=True)
    name = Column(String, nullable=False, unique=True)
    date = Column(DateTime, nullable=True)

    def __init__(self, lab_id=None, name=None, date=None):
        super().__init__()
        self.lab_id = lab_id
        self.name = name
        self.date = date

    def __repr__(self):
        return f"Run({self.lab_id!r}, {self.name!r}, {self.date!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Readset(BaseTable):
    """
    Readset:
        id integer [PK]
        sample_id integer [ref: > sample.id]
        run_id integer [ref: > run.id]
        name text
        lane text
        adapter1 text
        adapter2 text
        sequencing_type text
        quality_offset text
        alias blob
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "readset"

    sample_id = Column(Integer, ForeignKey("sample.id"))
    run_id = Column(Integer, ForeignKey("run.id"))
    name = Column(String, nullable=False, unique=True)
    lane = Column(String, nullable=True)
    adapter1 = Column(String, nullable=True)
    adapter2 = Column(String, nullable=True)
    sequencing_type = Column(String, nullable=True)
    quality_offset = Column(String, nullable=True)
    alias = Column(String, nullable=True)

    sample = relationship("Sample", backref="readset", lazy=False)
    run = relationship("Run", backref="readset", lazy=False)

    def __init__(self, name=None, lane=None, adapter1=None, adapter2=None, sequencing_type=None, quality_offset=None, alias=None):
        super().__init__()
        self.name = name
        self.lane = lane
        self.adapter1 = adapter1
        self.adapter2 = adapter2
        self.sequencing_type = sequencing_type
        self.quality_offset = quality_offset
        self.alias = alias

    def __repr__(self):
        return f"Readset({self.sample!r}, {self.run!r}, {self.name!r}, {self.lane!r}, {self.adapter1!r}, {self.adapter2!r}, {self.sequencing_type!r}, {self.quality_offset!r}, {self.alias!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Step(BaseTable):
    """
    Step:
        id integer [PK]
        sample_id integer [ref: > sample.id]
        readset_id integer [ref: > readset.id]
        name text
        status text
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "step"

    sample_id = Column(Integer, ForeignKey("sample.id"))
    readset_id = Column(Integer, ForeignKey("readset.id"))
    name = Column(String, nullable=False)
    status = Column(String, nullable=True)

    sample = relationship("Sample", backref="step", lazy=False)
    readset = relationship("Readset", backref="step", lazy=False)

    def __init__(self, name=None, status=None):
        super().__init__()
        self.name = name
        self.status = status

    def __repr__(self):
        return f"Step({self.sample!r}, {self.readset!r}, {self.name!r}, {self.status!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Job(BaseTable):
    """
    Job:
        id integer [PK]
        step_id integer [ref: > step.id]
        name text
        start timestamp
        stop timestamp
        status text
        type text
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "job"

    step_id = Column(Integer, ForeignKey("step.id"))
    name = Column(String, nullable=False)
    start = Column(DateTime, nullable=True)
    stop = Column(DateTime, nullable=True)
    status = Column(String, nullable=True)
    type = Column(String, nullable=True)

    step = relationship("Step", backref="job", lazy=False)

    def __init__(self, name=None, start=None, stop=None, status=None, type=None):
        super().__init__()
        self.name = name
        self.start = start
        self.stop = stop
        self.status = status
        self.type = type

    def __repr__(self):
        return f"Job({self.step!r}, {self.name!r}, {self.start!r}, {self.stop!r}, {self.status!r}, {self.type!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Metric(BaseTable):
    """
    Metric:
        id integer [PK]
        job_id integer [ref: > job.id]
        name text
        value text
        flag text //pass, warn, fail
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "metric"

    job_id = Column(Integer, ForeignKey("job.id"))
    name = Column(String, nullable=False)
    value = Column(String, nullable=True)
    flag = Column(String, nullable=True)

    job = relationship("Job", backref="metric", lazy=False)

    def __init__(self, name=None, value=None, flag=None):
        super().__init__()
        self.name = name
        self.value = value
        self.flag = flag

    def __repr__(self):
        return f"Metric({self.job!r}, {self.name!r}, {self.value!r}, {self.flag!r}, {self.deleted!r}, {self.extra_metadata!r})"

class File(BaseTable):
    """
    File:
        id integer [PK]
        job_id integer [ref: > job.id]
        path text
        type text
        description text
        creation timestamp
        deleted boolean
        extra_metadata json
    """
    __tablename__ = "file"

    job_id = Column(Integer, ForeignKey("job.id"))
    path = Column(String, nullable=True)
    type = Column(String, nullable=True)
    description = Column(String, nullable=True)
    creation = Column(DateTime, nullable=True)

    job = relationship("Job", backref="file", lazy=False)

    def __init__(self, path=None, type=None, description=None, creation=None):
        super().__init__()
        self.path = path
        self.type = type
        self.description = description
        self.creation = creation

    def __repr__(self):
        return f"File({self.job!r}, {self.path!r}, {self.type!r}, {self.description!r}, {self.creation!r}, {self.deleted!r}, {self.extra_metadata!r})"

def insert(engine, entry, **relations):
    """
    entry: Table object ex. Movie("Robocop", 1987)
    **relations: name = value ex. director = Director("Paul Verhoeven")
    """
    local_session = sessionmaker(autoflush=False, autocommit=False, bind=engine)
    # With this we get a session to do whatever we want to do
    session = local_session()

    # new_entry = entry
    for name, relation_entry in relations.items():
        setattr(entry, name, relation_entry)
        # print(getattr(relation_entry, "name"))
        # print(session.query(type(relation_entry)).populate_existing().all())
        # exit()
        # print(type(relation_entry).__name__)
        # print(session.query(type(relation_entry)).first())
        # if session.query(type(relation_entry)).first():
        #     setattr(entry, name, session.query(type(relation_entry)).first())
        # new_entry.name = relation_entry

    try:
        session.add(entry)
        session.commit()
    except Exception as error:
        print(f"Error: {error}")
        session.rollback()

# Cf. https://stackoverflow.com/questions/48722835/custom-type-hint-annotation
# T = typing.TypeVar('T')

# class json(typing.Generic[T]):
#     pass

# class timestamp(typing.Generic[T]):
#     pass
