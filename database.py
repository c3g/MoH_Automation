"""Module providing database tables and operations support."""
from sqlalchemy import Column, ForeignKey, Integer, Boolean, String, JSON, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relation, sessionmaker

Base = declarative_base()

class BaseTable(Base):
    """docstring for BaseTable"""

    id = Column(Integer, primary_key=True)
    deleted = Column(Boolean, default=False)
    extra_metadata = Column(JSON, nullable=True)

    def __init__(self, deleted=None, extra_metadata=None):
        self.deleted = deleted
        self.extra_metadata = extra_metadata

class Project(BaseTable):
    """docstring for Project"""
    __tablename__ = "project"

    name = Column(String, nullable=False)

    def __init__(self, name=None):
        super().__init__()
        self.name = name

    def __repr__(self):
        return f"Project({self.name!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Patient(BaseTable):
    """docstring for Patient"""
    __tablename__ = "patient"

    project_id = Column(Integer, ForeignKey("project.id"))
    name = Column(String, nullable=False, unique=True)
    alias = Column(String, nullable=True)
    cohort = Column(String, nullable=True)
    institution = Column(String, nullable=True)

    project = relation("Project", backref="patient", lazy=False)

    def __init__(self, name=None, alias=None, cohort=None, institution=None):
        super().__init__()
        self.name = name
        self.alias = alias
        self.cohort = cohort
        self.institution = institution

    def __repr__(self):
        return f"Patient({self.project!r}, {self.name!r}, {self.alias!r}, {self.cohort!r}, {self.institution!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Experiment(BaseTable):
    """docstring for experiment"""
    __tablename__ = "experiment"

    run_id = Column(Integer, ForeignKey("run.id"))
    project_id = Column(Integer, ForeignKey("project.id"))
    sequencing_technology = Column(String, nullable=True)

    run = relation("Run", backref="experiment", lazy=False)
    project = relation("Project", backref="experiment", lazy=False)

    def __init__(self, sequencing_technology=None):
        super().__init__()
        self.sequencing_technology = sequencing_technology

    def __repr__(self):
        return f"Experiment({self.run!r}, {self.project!r}, {self.sequencing_technology!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Sample(BaseTable):
    """docstring for Sample"""
    __tablename__ = "sample"

    patient_id = Column(Integer, ForeignKey("patient.id"))
    experiment_id = Column(Integer, ForeignKey("experiment.id"))
    name = Column(String, nullable=False, unique=True)
    tumour = Column(Boolean, default=False)
    alias = Column(String, nullable=True)

    patient = relation("Patient", backref="sample", lazy=False)
    experiment = relation("Patient", backref="sample", lazy=False)

    def __init__(self, name=None, tumour=None, alias=None):
        super().__init__()
        self.name = name
        self.tumour = tumour
        self.alias = alias

    def __repr__(self):
        return f"Sample({self.patient!r}, {self.experiment!r}, {self.name!r}, {self.tumour!r}, {self.alias!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Run(BaseTable):
    """docstring for Patient"""
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
    """docstring for Readset"""
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

    sample = relation("Sample", backref="readset", lazy=False)
    run = relation("Run", backref="readset", lazy=False)

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
    """docstring for Step"""
    __tablename__ = "step"

    sample_id = Column(Integer, ForeignKey("sample.id"))
    readset_id = Column(Integer, ForeignKey("readset.id"))
    name = Column(String, nullable=False)
    status = Column(String, nullable=True)

    sample = relation("Sample", backref="step", lazy=False)
    readset = relation("Readset", backref="step", lazy=False)

    def __init__(self, name=None, status=None):
        super().__init__()
        self.name = name
        self.status = status

    def __repr__(self):
        return f"Step({self.sample!r}, {self.readset!r}, {self.name!r}, {self.status!r}, {self.deleted!r}, {self.extra_metadata!r})"

class Job(BaseTable):
    """docstring for Job"""
    __tablename__ = "job"

    step_id = Column(Integer, ForeignKey("step.id"))
    name = Column(String, nullable=False)
    start = Column(DateTime, nullable=True)
    stop = Column(DateTime, nullable=True)
    status = Column(String, nullable=True)
    type = Column(String, nullable=True)

    step = relation("Step", backref="job", lazy=False)

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
    """docstring for Metric"""
    __tablename__ = "metric"

    job_id = Column(Integer, ForeignKey("job.id"))
    name = Column(String, nullable=False)
    value = Column(String, nullable=True)
    flag = Column(String, nullable=True)

    job = relation("Job", backref="metric", lazy=False)

    def __init__(self, name=None, value=None, flag=None):
        super().__init__()
        self.name = name
        self.value = value
        self.flag = flag

    def __repr__(self):
        return f"Metric({self.job!r}, {self.name!r}, {self.value!r}, {self.flag!r}, {self.deleted!r}, {self.extra_metadata!r})"

class File(BaseTable):
    """docstring for File"""
    __tablename__ = "file"

    job_id = Column(Integer, ForeignKey("job.id"))
    path = Column(String, nullable=True)
    type = Column(String, nullable=True)
    description = Column(String, nullable=True)
    creation = Column(DateTime, nullable=True)

    job = relation("Job", backref="file", lazy=False)

    def __init__(self, path=None, type=None, description=None, creation=None):
        super().__init__()
        self.path = path
        self.type = type
        self.description = description
        self.creation = creation

    def __repr__(self):
        return f"File({self.job!r}, {self.path!r}, {self.type!r}, {self.description!r}, {self.creation!r}, {self.deleted!r}, {self.extra_metadata!r})"

# Cf. https://stackoverflow.com/questions/48722835/custom-type-hint-annotation
# T = typing.TypeVar('T')

# class json(typing.Generic[T]):
#     pass

# class timestamp(typing.Generic[T]):
#     pass
