from attrs import define, field
from solid2.core.object_base import BareOpenSCADObject
from solid2 import cube, cylinder, sphere
from math import sin, cos, pi
import numpy as np
from collections import namedtuple
from typing import Annotated
import numpy.typing as npt


DTR = 360/2/pi
PRIM = "solid2.core.builtins.openscad_primitives."
BOSL = "solid2.extensions.bosl2."


# # Example:
# from solidbox import Bbox
# from solid2 import cube
# some_scad = cube([1,1,1]).rotate([-20, 20, 10]).scale([2,1,1]) + cube([2,1,1]).up(10)
# bbox_frame_scad = Bbox.from_scad(some_scad).as_frame

Point3D = namedtuple("Point3D", ["x", "y", "z"])
Mat = Annotated[npt.NDArray[np.float64], (3, 3)]
MAT_ID: Mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


@define
class Bbox:
    x_min: float = field(default=0, converter=float)
    x_max: float = field(default=0, converter=float)
    y_min: float = field(default=0, converter=float)
    y_max: float = field(default=0, converter=float)
    z_min: float = field(default=0, converter=float)
    z_max: float = field(default=0, converter=float)
    
    depth: float = field(init=False)
    width: float = field(init=False)
    height: float = field(init=False)
    center: Point3D = field(init=False)

    def __attrs_post_init__(self):
        self.depth = self.x_max - self.x_min
        self.width = self.y_max - self.y_min
        self.height = self.z_max - self.z_min
        self.center = Point3D((self.x_max + self.x_min)/2, (self.y_max + self.y_min)/2, (self.z_max + self.z_min)/2)
    
    @property
    def x(self):
        return self.depth

    @property
    def y(self):
        return self.width

    @property
    def z(self):
        return self.height

    @property
    def size(self):
        return self.x, self.y, self.z

    @property
    def center_n(self):
        return -self.center.x, -self.center.y, -self.center.z
    
    @property
    def is_null(self) -> bool:
        return all(v == 0 for v in (self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max))

    @property
    def as_tuple(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        return ((self.x_min, self.x_max), (self.y_min, self.y_max), (self.z_min, self.z_max))
    
    @property
    def points(
        self,
    ) -> tuple[
        Point3D,
        Point3D,
        Point3D,
        Point3D,
        Point3D,
        Point3D,
        Point3D,
        Point3D,
    ]:
        points = []
        for x in self.as_tuple[0]:
            for y in self.as_tuple[1]:
                for z in self.as_tuple[2]:
                    points.append((x, y, z))
        return tuple(points)

    @property
    def as_cube(self) -> BareOpenSCADObject:
        return cube(self.size, center=True).translate(self.center)

    @property
    def as_frame(self) -> BareOpenSCADObject:
        R = min(self.x, self.y, self.z) / 50
        sides = (
            cylinder(r=R, h=self.z).translate([self.x_min, self.y_min, self.z_min]) +
            cylinder(r=R, h=self.z).translate([self.x_min, self.y_max, self.z_min]) +
            cylinder(r=R, h=self.z).translate([self.x_max, self.y_min, self.z_min]) +
            cylinder(r=R, h=self.z).translate([self.x_max, self.y_max, self.z_min]) 
        )
        bottom = (
            cylinder(r=R, h=self.y).rotateX(-90).translate([self.x_min, self.y_min, self.z_min]) + 
            cylinder(r=R, h=self.y).rotateX(-90).translate([self.x_max, self.y_min, self.z_min]) + 
            cylinder(r=R, h=self.x).rotateX(-90).rotateZ(-90).translate([self.x_min, self.y_min, self.z_min]) + 
            cylinder(r=R, h=self.x).rotateX(-90).rotateZ(-90).translate([self.x_min, self.y_max, self.z_min])
        )
        d = sides + bottom + bottom.up(self.z)
        for x in self.as_tuple[0]:
            for y in self.as_tuple[1]:
                for z in self.as_tuple[2]:
                    d += sphere(r=R).translate(x, y, z)
        return d
    
    @classmethod
    def from_tuple(cls, tup: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]) -> "Bbox":
        xt, yt, zt = tuple(tup)
        xn, xx = xt
        yn, yx = yt
        zn, zx = zt
        return cls(xn, xx, yn, yx, zn, zx)
    
    @classmethod
    def from_points(cls, *points: Point3D) -> "Bbox":
        xes = [p.x for p in points]
        yes = [p.y for p in points]
        zes = [p.z for p in points]
        return cls(min(xes), max(xes), min(yes), max(yes), min(zes), max(zes))
    
    @classmethod
    def from_scad(cls, d: BareOpenSCADObject, expedite: Mat = MAT_ID, quiet: bool = False) -> "Bbox":
        name = f"{type(d).__module__}.{type(d).__name__}"
        params = d._params

        # Translate some operators
        if name == BOSL + "transforms.up":
            assert "z" in params
            name = PRIM + "translate"
            params = {"v": [0, 0, params["z"]]}
        elif name == BOSL + "transforms.down":
            assert "z" in params
            name = PRIM + "translate"
            params = {"v": [0, 0, -params["z"]]}
        elif name == BOSL + "transforms.right":
            assert "x" in params
            name = PRIM + "translate"
            params = {"v": [params["x"], 0, 0]}
        elif name == BOSL + "transforms.left":
            assert "x" in params
            name = PRIM + "translate"
            params = {"v": [-params["x"], 0, 0]}

        # Handle objects
        if name in [
            PRIM + "cube",
            BOSL + "shapes3d.cuboid",
        ]:
            # Bosl cuboids are centered by default (or always?), normal cubes are not
            if name == PRIM + "cube":
                c = params["center"] or False
            elif name == BOSL + "shapes3d.cuboid":
                c = params.get("center", True)
            else:
                assert False
            assert isinstance(c, bool)
            
            x, y, z = params["size"]
            obj = cls(-c/2*x, x-c/2*x, -c/2*y, y-c/2*y, -c/2*z, z-c/2*z)
            return obj._apply_matrix(expedite)

        elif name == PRIM + "translate":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            trans = expedite @ np.array(params["v"])
            return cls.from_tuple(
                (mn + tdim, mx + tdim)
                for tdim, (mn, mx) in zip(trans, cls.from_scad(d._children[0], expedite, quiet=quiet).as_tuple)  # pyright: ignore[reportArgumentType]
            )
        
        elif name == PRIM + "rotate":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            _a: tuple[float, float, float] = params["a"]  # pyright: ignore[reportAssignmentType]

            # This would return the rotated bbox of the union of children. It returns a valid upper bound for the box, but 
            # is needlessly conservative.
            # return Bbox.from_scad(d._children[0], quiet=quiet).__rotate(_a)

            # Instead, try to do the rotation as early as possible
            return Bbox.from_scad(d._children[0], expedite @ rot.xyz(_a), quiet=quiet)

        elif name == PRIM + "scale":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            _a: tuple[float, float, float] = params["v"]  # pyright: ignore[reportAssignmentType]
            return Bbox.from_scad(d._children[0], expedite @ scl.xyz(_a), quiet=quiet)

        elif name in [
            PRIM + "union",
            BOSL + "regions.union",
            BOSL + "color.hsv",
        ]:
            if boxes := [box for box in [cls.from_scad(child, expedite, quiet=quiet) for child in d._children] if not box.is_null]:
                return cls._union(*boxes)
            else:
                # No children or all were null, this branch does not contribute to Bbox.
                return cls()

        elif name in [
            PRIM + "difference",
            BOSL + "regions.difference",
        ]:
            # NOTE: This is an upper bound. It should be possible to compute the difference of all bboxes.
            return cls.from_scad(d._children[0], expedite, quiet=quiet)

        elif name in [
            PRIM + "cylinder",
            BOSL + "shapes3d.cylinder",
        ]:
            h = params["h"]
            rmax = params["r"] or max(params[rr] or 0 for rr in ("r1", "r2"))
            assert rmax != 0
            # TODO: implement diameter (d, d1, d2)
            c = int(params["center"] or 0)  # pyright: ignore[reportArgumentType]
            obj = cls(-rmax, rmax, -rmax, rmax, -c/2*h, h*(1-c) + c/2*h)  # pyright: ignore[reportArgumentType, reportOperatorIssue]
            return obj._apply_matrix(expedite)

        elif name in [
            PRIM + "sphere",
            BOSL + "shapes3d.sphere",
        ]:
            r = params["r"]
            obj = cls(-r, r, -r, r, -r, r)  # pyright: ignore[reportArgumentType, reportOperatorIssue]
            return obj._apply_matrix(expedite)

        elif name in [
            BOSL + "regions.intersection"
        ]:
            if boxes := [box for box in [cls.from_scad(child, expedite, quiet=quiet) for child in d._children] if not box.is_null]:
                return cls._intersection(*boxes)
            else:
                # No children or all were null, this branch does not contribute to Bbox.
                return cls()

        else:
            if not quiet:
                print(f"Warning: Unsupported object; '{name}'. Ignoring.")
            # Return 'special' Bbox for which is_null will be true (all-zero).
            # Such a Bbox will be ignored by further Bbox computations.
            # IDEA: perhaps I could have the user add a bbox annotation to some objects they know aren't implemented/able. Could be as simple as just setting a 3x2-tuple to some new attribute, e.g. r=rack(...); r.bbox_hint=(...). This hint could then be used later when computing the bbox if the object type is not matched.
            return cls()

    
    def _apply_matrix(self, mat: Mat) -> "Bbox":
        return self.from_points(*(Point3D(*(mat @ p)) for p in self.points))

    @classmethod
    def __union_intersection(cls, minmaxmaxmin, *boxes: "Bbox") -> "Bbox":
        return cls.from_tuple(
            [
                minmaxmaxmin[inx](box.as_tuple[idim][inx] for box in boxes)
                for inx in (0, 1)
            ]
            for idim in (0, 1, 2)  # pyright: ignore[reportArgumentType]
        )

    @classmethod
    def _union(cls, *boxes: "Bbox") -> "Bbox":
        return cls.__union_intersection((min, max), *boxes)

    @classmethod
    def _intersection(cls, *boxes: "Bbox") -> "Bbox":
        return cls.__union_intersection((max, min), *boxes)


class mat:
    @classmethod
    def xyz(cls, vector: Point3D | tuple[float, float, float]) -> Mat:
        bx, by, bz = vector
        return cls.z(bz) @ cls.y(by) @ cls.x(bx)  # pyright: ignore[reportAttributeAccessIssue]


class rot(mat):
    @staticmethod
    def x(theta: float) -> Mat:
        theta /= DTR
        return np.array([
            [1, 0, 0],
            [0, cos(theta), -sin(theta)],
            [0, sin(theta), cos(theta)],
        ])

    @staticmethod
    def y(theta: float) -> Mat:
        theta /= DTR
        return np.array([
            [cos(theta), 0, sin(theta)],
            [0, 1, 0],
            [-sin(theta), 0, cos(theta)],
        ])

    @staticmethod
    def z(theta: float) -> Mat:
        theta /= DTR
        return np.array([
            [cos(theta), -sin(theta), 0],
            [sin(theta), cos(theta), 0],
            [0, 0, 1],
        ])
    

class scl(mat):
    @staticmethod
    def x(s: float) -> Mat:
        return np.array([
            [s, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ])

    @staticmethod
    def y(s: float) -> Mat:
        return np.array([
            [1, 0, 0],
            [0, s, 0],
            [0, 0, 1],
        ])

    @staticmethod
    def z(s: float) -> Mat:
        return np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, s],
        ])
