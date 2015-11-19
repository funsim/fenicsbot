from dolfin import *
from mshr import *
import dolfin

def UnitSquareMesh(res=10):
    return dolfin.UnitSquareMesh(10, 10)

def DolfinMesh():
    return Mesh("dolfin.xml.gz")

def GearMesh():
    return Mesh("gear.xml.gz")

def HoleMesh():
    domain = Rectangle(dolfin.Point(0., 0.), dolfin.Point(1., 1.)) - \
             Circle(dolfin.Point(0.0, 0.0), .35)
    return generate_mesh(domain, 32)


if __name__ == "__main__":
    plot(DolfinMesh())
    plot(UnitSquareMesh())
    interactive()

