from abc import ABCMeta, abstractmethod
from attrs import define, field
from solid2.core.object_base import BareOpenSCADObject
from .main import Bbox


@define(repr=False)
class _Util(metaclass=ABCMeta):
    _scad_object: BareOpenSCADObject
    _quiet: bool = True
    bbox: Bbox = field(init=False)
    
    def __attrs_post_init__(self):
        self.bbox = Bbox.from_scad(self._scad_object, quiet=self._quiet)

    @property
    @abstractmethod
    def x(self) -> float:
        ...

    @property
    @abstractmethod
    def y(self) -> float:
        ...

    @property
    @abstractmethod
    def z(self) -> float:
        ...

    def __repr__(self, attrs=("x", "y", "z")):
        dims = [f"{p} = {getattr(self, p)}" for p in attrs]
        return type(self).__name__ + "\n  ".join(("", *dims))


@define(repr=False)
class Size(_Util):
    @property
    def x(self):
        return self.bbox.depth

    @property
    def y(self):
        return self.bbox.width

    @property
    def z(self):
        return self.bbox.height


@define(repr=False)
class Mid(_Util):
    @property
    def x(self):
        return self.bbox.center.x

    @property
    def y(self):
        return self.bbox.center.y

    @property
    def z(self):
        return self.bbox.center.z


@define(repr=False)
class Max(_Util):
    @property
    def x(self):
        return self.bbox.x_max

    @property
    def y(self):
        return self.bbox.y_max

    @property
    def z(self):
        return self.bbox.z_max


@define(repr=False)
class Min(_Util):
    @property
    def x(self):
        return self.bbox.x_min

    @property
    def y(self):
        return self.bbox.y_min

    @property
    def z(self):
        return self.bbox.z_min
