from dolfin import *

def DolfinMesh():
    return Mesh("dolfin.xml.gz")

def GearMesh():
    return Mesh("gear.xml.gz")

if __name__ == "__main__":
    plot(DolfinMesh())
    interactive()

