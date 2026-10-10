# Bracket-style multi-feature script: two bolt holes, a fillet attempt,
# element-type override, two steps (tip pull + self weight + edge
# pressure). Exercises the unsupported/warning paths. SI units.
from abaqus import *
from abaqusConstants import *
import mesh

s = mdb.models['Model-1'].ConstrainedSketch(name='profile', sheetSize=1.0)
s.rectangle(point1=(0.0, 0.0), point2=(0.5, 0.3))
s.Line(point1=(0.0, 0.3), point2=(0.5, 0.3))
s.CircleByCenterRadius(center=(0.1, 0.15), radius=0.03)
s.CircleByCenterRadius(center=(0.4, 0.15), radius=0.03)
s.Fillet(radius=0.02, vertex1=(0.0, 0.0), vertex2=(0.5, 0.0))
p = mdb.models['Model-1'].Part(name='Bracket', dimensionality=TWO_D_PLANAR,
    type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)

mdb.models['Model-1'].Material(name='Alu')
mdb.models['Model-1'].materials['Alu'].Elastic(table=((70000000000.0, 0.33), ))
mdb.models['Model-1'].materials['Alu'].Density(table=((2700.0,), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Sec', material='Alu',
    thickness=0.005)
p.SectionAssignment(region=p.Set(faces=p.faces, name='All'),
    sectionName='Sec')

elemType = mesh.ElemType(elemCode=CPS8, elemLibrary=STANDARD)
p.setElementType(regions=(p.faces,), elemTypes=(elemType,))

p.seedPart(size=0.02)
p.generateMesh()

a = mdb.models['Model-1'].rootAssembly
a.Instance(name='Bracket-1', part=p, dependent=ON)

mdb.models['Model-1'].StaticStep(name='Load', previous='Initial')

fix = a.instances['Bracket-1'].edges.findAt(((0.0, 0.15, 0.0), ))
mdb.models['Model-1'].EncastreBC(name='Fix', createStepName='Load', region=fix)

pull = a.instances['Bracket-1'].vertices.findAt(((0.5, 0.3, 0.0), ))
mdb.models['Model-1'].ConcentratedForce(name='Pull', createStepName='Load',
    region=pull, cf1=500.0, cf2=-200.0)

mdb.models['Model-1'].StaticStep(name='Gravity', previous='Load')
mdb.models['Model-1'].GravityLoad(name='SelfWeight', createStepName='Gravity',
    gravitation=(0.0, -9.81, 0.0))

topEdge = a.instances['Bracket-1'].edges.findAt(((0.25, 0.3, 0.0), ))
mdb.models['Model-1'].Pressure(name='Top', createStepName='Gravity',
    region=topEdge, magnitude=250000.0)

mdb.save()
