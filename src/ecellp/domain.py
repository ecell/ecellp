__program__ = 'domain'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu <kaizu@riken.jp>'
__copyright__ = ''
__license__ = ''


class Domain(object):

    def __init__(self, name, **attrs):
        self.__name = name
        self.__attributes = {}

        self.set_attributes(**attrs)

    def set_attributes(self, **attrs):
        if "name" in attrs.keys():
            raise RuntimeError, 'Do not overload "name"'
        self.__attributes.update(attrs)

    def get_attribute(self, key):
        return self.__attributes.get(key)

    def attributes(self):
        return self.__attributes

    def serial(self):
        return self.name()

    def name(self):
        return self.__name

    def match(self, dom, parent):
        return dom.name() == self.name()

    def as_Domain(self):
        return Domain(self.name(), **self.attributes())

    def __repr__(self):
        mod = self.__class__.__module__
        cls = self.__class__.__name__
        # mem = '0x' + hex(id(self))[2:].zfill(8).upper()
        return '<{0}.{1} instance with serial "{2}">'.format(mod, cls, self.serial())

class FilteringDomain(Domain):

    def __init__(self, name = "", **attrs):
        Domain.__init__(self, name, **attrs)

    def match(self, dom, parent):
        if self.name() is not "" and self.name() != dom.name():
            return False

        for key, value in self.attributes().items():
            if dom.get_attribute(key) != value:
                return False
        else:
            return True

class Condition(object):

    def __init__(self, left, right, includes = [], excludes = []):
        if len(includes) == 0 and len(excludes) == 0:
            raise RuntimeError, ""

        self.__bounds = (left, right)
        self.__includes = includes
        self.__excludes = excludes

    def match(self, dom, parent):
        start = dom.get_attribute("start")
        stride = dom.get_attribute("stride")
        region = (
            start + self.__bounds[0] * stride,
            start + self.__bounds[1] * stride,
            stride)
        neighbors = parent.query_domains_by_region(*region)

        # do more precise filtering before checking

        for elem in self.__includes:
            if not any([elem.match(neighbor, parent) for neighbor in neighbors]):
                return False
        for elem in self.__excludes:
            if any([elem.match(neighbor, parent) for neighbor in neighbors]):
                return False
        return True

class ConditionalDomain(Domain):

    def __init__(self, name, **attrs):
        Domain.__init__(self, name, **attrs)

        self.__conditions = []

    def add_condition(self, condition):
        self.__conditions.append(condition)

    def find(self, parent):
        retval = []
        for dom in parent.list_domains(self.as_Domain()):
            if self.match(dom, parent):
                retval.append(dom)
        return retval

    def match(self, dom, parent):
        if not self.as_Domain().match(dom, parent):
            return False

        for condition in self.__conditions:
            if not condition.match(dom, parent):
                return False
        else:
            return True
