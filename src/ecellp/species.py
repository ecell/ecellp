from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

class BaseDec(Base):
    __tablename__ = 'base'
    __table_args__ = {'sqlite_autoincrement': True}

    dectype = Column(String, nullable=False)
    __mapper_args__ = {'polymorphic_on': dectype}

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    strand = Column(Integer, nullable=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    feature = Column(String, nullable=False)
    sequence = Column(String, nullable=False)

    def __init__(self, name, strand, start, end, feature, sequence):
        self.name = name
        self.strand = int(strand)
        self.start = int(start)
        self.end = int(end)
        self.feature = feature
        self.sequence = sequence

    def keys(self):
        return (
            "name", "strand", "start", "end", "feature", "sequence")

    def values(self):
        return (
            self.name, self.strand, self.start, self.end,
            self.feature, self.sequence)

    def as_dict(self):
        retval = {}
        for key, value in zip(self.keys(), self.values()):
            retval[key] = value
        return retval

    def __repr__(self):
        return "<%s.%s: %s>" % (
            self.__class__.__module__, self.__class__.__name__,
            str(self.as_dict()))

class CDSDec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'cds'}

    def __init__(self, *args):
        super(CDSDec, self).__init__(*args)

class tRNADec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'trna'}

    def __init__(self, *args):
        super(tRNADec, self).__init__(*args)

class rRNADec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'rrna'}

    def __init__(self, *args):
        super(rRNADec, self).__init__(*args)

class PromoterDec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'promoter'}

    def __init__(self, *args):
        super(PromoterDec, self).__init__(*args)

class TerminatorDec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'terminator'}

    def __init__(self, *args):
        super(TerminatorDec, self).__init__(*args)

class OperonDec(BaseDec):
    __mapper_args__ = {'polymorphic_identity': 'operon'}

    def __init__(self, *args):
        super(TerminatorDec, self).__init__(*args)
