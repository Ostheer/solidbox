from attrs import define, field
from solid2.core.object_base.object_base_impl import BareOpenSCADObject
from solid2 import cube, cylinder, sphere
from math import sin, cos, pi
from typing import Callable, Mapping, Sequence
from scipy.spatial.transform import Rotation as R

DTR = 360/2/pi
PRIM = "solid2.core.builtins.openscad_primitives."
BOSL = "solid2.extensions.bosl2."


# # Example:
# from solidbox import Bbox
# from solid2 import cube
# some_scad = cube([1,1,1]).rotate([-20, 20, 10]).scale([2,1,1]) + cube([2,1,1]).up(10)
# bbox_frame_scad = Bbox.from_scad(some_scad).as_frame

def apply(obj: "Bbox", *operations: "Operation") -> "Bbox":
    for op in operations:
        obj = op.callback(obj, **op.kwargs)
    return obj


@define
class Bbox:
    x_min: float = 0
    x_max: float = 0
    y_min: float = 0
    y_max: float = 0
    z_min: float = 0
    z_max: float = 0
    
    depth: float = field(init=False)
    width: float = field(init=False)
    height: float = field(init=False)
    center: tuple[float, float, float] = field(init=False)

    def __attrs_post_init__(self):
        self.depth = self.x_max - self.x_min
        self.width = self.y_max - self.y_min
        self.height = self.z_max - self.z_min
        self.center = (self.x_max + self.x_min)/2, (self.y_max + self.y_min)/2, (self.z_max + self.z_min)/2
    
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
        return -self.center[0], -self.center[1], -self.center[2]
    
    @property
    def is_null(self) -> bool:
        return all(v == 0 for v in (self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max))

    @property
    def as_tuple(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        return ((self.x_min, self.x_max), (self.y_min, self.y_max), (self.z_min, self.z_max))

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
    def from_tuple(cls, tup) -> "Bbox":
        tup = tuple(tup)
        xn, xx = tup[0]
        yn, yx = tup[1]
        zn, zx = tup[2]
        return cls(xn, xx, yn, yx, zn, zx)
    
    @classmethod
    def from_scad(cls, d: BareOpenSCADObject, *expedite: "Operation") -> "Bbox":
        def apply_expedited(d):
            return apply(d, *expedite)
        
        name = f"{type(d).__module__}.{type(d).__name__}"
        params = d._params

        expedite = merge_operations(*expedite)  # pyright: ignore[reportAssignmentType]

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
            return apply_expedited(obj)

        elif name == PRIM + "translate":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            trans: tuple[float, float, float] = params["v"]  # pyright: ignore[reportAssignmentType]
            for op in expedite:
                # The constituent objects get rotated below, but here we also rotate and scale the translation itself.
                if isinstance(op, Rotation):
                    trans = rotate_point(trans, op.vector)
                elif isinstance(op, Scaling):
                    trans = scale_point(trans, op.vector)
            return cls.from_tuple(
                (mn + tdim, mx + tdim)
                for tdim, (mn, mx) in zip(trans, cls.from_scad(d._children[0], *expedite).as_tuple)
            )
        
        elif name == PRIM + "rotate":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            _a: tuple[float, float, float] = params["a"]  # pyright: ignore[reportAssignmentType]

            # This would return the rotated bbox of the union of children. It returns a valid upper bound for the box, but 
            # is needlessly conservative.
            # return Bbox.from_scad(d._children[0]).__rotate(_a)

            # Instead, try to do the rotation as early as possible
            return Bbox.from_scad(d._children[0], Rotation(callback=cls._rotate, kwargs={"vector": _a}), *expedite)

        elif name == PRIM + "scale":
            assert len(d._children) == 1  # Is there always either a primitive or a union here? Else must unionize the children.
            _a: tuple[float, float, float] = params["v"]  # pyright: ignore[reportAssignmentType]
            return Bbox.from_scad(d._children[0], Scaling(callback=cls._scale, kwargs={"vector": _a}), *expedite)

        elif name in [
            PRIM + "union",
            BOSL + "regions.union",
            BOSL + "color.hsv",
        ]:
            if boxes := [box for box in [cls.from_scad(child, *expedite) for child in d._children] if not box.is_null]:
                return cls._union(*boxes)
            else:
                # No children or all were null, this branch does not contribute to Bbox.
                return cls()

        elif name in [
            PRIM + "difference",
            BOSL + "regions.difference",
        ]:
            # NOTE: This is an upper bound. It should be possible to compute the difference of all bboxes.
            return cls.from_scad(d._children[0], *expedite)

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
            return apply_expedited(obj)

        elif name in [
            BOSL + "regions.intersection"
        ]:
            if boxes := [box for box in [cls.from_scad(child, *expedite) for child in d._children] if not box.is_null]:
                return cls._intersection(*boxes)
            else:
                # No children or all were null, this branch does not contribute to Bbox.
                return cls()

        else:
            print(f"Warning: Unsupported object; '{name}'. Ignoring.")
            # Return 'special' Bbox for which is_null will be true (all-zero).
            # Such a Bbox will be ignored by further Bbox computations.
            # IDEA: perhaps I could have the user add a bbox annotation to some objects they know aren't implemented/able. Could be as simple as just setting a 3x2-tuple to some new attribute, e.g. r=rack(...); r.bbox_hint=(...). This hint could then be used later when computing the bbox if the object type is not matched.
            return cls()


    def _rotate(
        self,
        vector: tuple[float, float, float],
        ):
        rotated_bbox_points_x = []
        rotated_bbox_points_y = []
        rotated_bbox_points_z = []
        for xo in self.as_tuple[0]:
            for yo in self.as_tuple[1]:
                for zo in self.as_tuple[2]:
                    x, y, z = rotate_point((xo, yo, zo), vector)
                    rotated_bbox_points_x.append(x)
                    rotated_bbox_points_y.append(y)
                    rotated_bbox_points_z.append(z)

        return type(self)(
            min(rotated_bbox_points_x),
            max(rotated_bbox_points_x),
            min(rotated_bbox_points_y),
            max(rotated_bbox_points_y),
            min(rotated_bbox_points_z),
            max(rotated_bbox_points_z),
        )

    def _scale(
        self,
        vector: tuple[float, float, float],
        ):
        xmin, ymin, zmin = scale_point((self.x_min, self.y_min, self.z_min), vector)
        xmax, ymax, zmax = scale_point((self.x_max, self.y_max, self.z_max), vector)
        return type(self)(
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax,
        )

    @classmethod
    def __union_intersection(cls, minmaxmaxmin, *boxes: "Bbox") -> "Bbox":
        return cls.from_tuple(
            [
                minmaxmaxmin[inx](box.as_tuple[idim][inx] for box in boxes)
                for inx in (0, 1)
            ]
            for idim in (0, 1, 2)
        )

    @classmethod
    def _union(cls, *boxes: "Bbox") -> "Bbox":
        return cls.__union_intersection((min, max), *boxes)

    @classmethod
    def _intersection(cls, *boxes: "Bbox") -> "Bbox":
        return cls.__union_intersection((max, min), *boxes)


def rotate_point(
    point: tuple[float, float, float], 
    rotation: tuple[float, float, float]
    ) -> tuple[float, float, float]:
    x, y, z = point
    bx, by, bz = rotation

    # Rotate around X-axis
    # Treat Y-component
    y1 = y * cos(bx / DTR)
    z1 = y * sin(bx / DTR)
    # Treat Z-component
    z2 = z * cos(bx / DTR)
    y2 = z * -sin(bx / DTR)
    # Sum
    y = y1 + y2
    z = z1 + z2

    # Rotate around Y-axis
    # Treat Z-component
    z1 = z * cos(by / DTR)
    x1 = z * sin(by / DTR)
    # Treat X-component
    x2 = x * cos(by / DTR)
    z2 = x * -sin(by / DTR)
    # Sum
    x = x1 + x2
    z = z1 + z2

    # Rotate around Z-axis
    # Treat X-component
    x1 = x * cos(bz / DTR)
    y1 = x * sin(bz / DTR)
    # Treat Y-component
    y2 = y * cos(bz / DTR)
    x2 = y * -sin(bz / DTR)
    # Sum
    x = x1 + x2
    y = y1 + y2

    return x, y, z


def scale_point(
    point: tuple[float, float, float], 
    scaling: tuple[float, float, float]
    ) -> tuple[float, float, float]:
    x, y, z = point
    bx, by, bz = scaling

    return bx*x, by*y, bz*z


@define
class Operation:
    kwargs: Mapping
    callback: Callable

    @property
    def vector(self) -> tuple[float, float, float]:
        return self.kwargs["vector"]

class Rotation(Operation):
    pass

class Scaling(Operation):
    pass


def merge_operations(*expedite: Operation) -> Sequence[Operation]:
        if not expedite:
            return ()
        
        new = []
        rot = None
        scl = None
        for op in expedite[::-1]:
            if isinstance(op, Rotation):
                # We have a rotation (maybe again)
                _rot = R.from_euler('xyz', op.vector, degrees=True)
                if rot is None:
                    rot = _rot
                else:
                    rot *= _rot
            
            elif rot is not None:
                # We had rotations, but now it's something else. Put hot rotations in the list
                _totrot = tuple([float(v) for v in rot.as_euler('xyz', degrees=True)])
                new.append(Rotation(callback=Bbox._rotate, kwargs={"vector": _totrot}),)
                rot = None

            if isinstance(op, Scaling):
                _scl = op.vector
                if scl is None:
                    scl = _scl
                else:
                    scl = [r1*r2 for r1, r2 in zip(scl, _scl)]
            
            elif scl is not None:
                new.append(Scaling(callback=Bbox._scale, kwargs={"vector": scl}),)
                scl = None

        # The last one still has to be added
        if rot is not None:
            _totrot = tuple([float(v) for v in rot.as_euler('xyz', degrees=True)])
            new.append(Rotation(callback=Bbox._rotate, kwargs={"vector": _totrot}),)
        elif scl is not None:
            new.append(Scaling(callback=Bbox._scale, kwargs={"vector": scl}),)
        
        return new[::-1]
