from solid2.core.object_base import BareOpenSCADObject
from .main import Bbox
from attrs import define, field


@define
class _Util:
    _scad_object: BareOpenSCADObject
    _quiet: bool = True
    bbox: Bbox = field(init=False)
    
    def __attrs_post_init__(self):
        self.bbox = Bbox.from_scad(self._scad_object, quiet=self._quiet)


@define
class Size(_Util):
    @property
    def x(self) -> float:
        return self.bbox.depth

    @property
    def y(self) -> float:
        return self.bbox.width

    @property
    def z(self) -> float:
        return self.bbox.height


@define
class Mid(_Util):
    @property
    def x(self) -> float:
        return self.bbox.center.x

    @property
    def y(self) -> float:
        return self.bbox.center.y

    @property
    def z(self) -> float:
        return self.bbox.center.z


@define
class Max(_Util):
    @property
    def x(self) -> float:
        return self.bbox.x_max

    @property
    def y(self) -> float:
        return self.bbox.y_max

    @property
    def z(self) -> float:
        return self.bbox.z_max


@define
class Min(_Util):
    @property
    def x(self) -> float:
        return self.bbox.x_min

    @property
    def y(self) -> float:
        return self.bbox.y_min

    @property
    def z(self) -> float:
        return self.bbox.z_min
