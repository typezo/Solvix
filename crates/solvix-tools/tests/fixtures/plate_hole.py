# Square plate with a central circular hole under uniaxial edge
# pressure — the Kirsch benchmark (far-field 1 MPa, hole-edge peak
# 3 MPa). SI units (m, Pa, N).
from abaqus import *
from abaqusConstants import *
import mesh

s = mdb.models['Model-1'].ConstrainedSketch(name='profile', sheetSize=4.0)
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-1.0, -1.0), point2=(1.0, 1.0))
s.CircleByCenterRadius(center=(0.0, 0.0), radius=0.2)
p = mdb.models['Model-1'].Part(name='Plate', dimensionality=TWO_D_PLANAR,
    type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)
s.unsetPrimaryObject()

mdb.models['Model-1'].Material(name='Steel')
mdb.models['Model-1'].materials['Steel'].Elastic(table=((210000000000.0, 0.3), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section', material='Steel',
    thickness=0.01)
p.SectionAssignment(region=p.Set(faces=p.faces, name='AllFaces'),
    sectionName='Section')

p.seedPart(size=0.05)
p.generateMesh()

a = mdb.models['Model-1'].rootAssembly
a.Instance(name='Plate-1', part=p, dependent=ON)

mdb.models['Model-1'].StaticStep(name='Load', previous='Initial')

leftEdge = a.instances['Plate-1'].edges.findAt(((-1.0, 0.0, 0.0), ))
mdb.models['Model-1'].DisplacementBC(name='FixX', createStepName='Load',
    region=leftEdge, u1=0.0)

corner = a.instances['Plate-1'].vertices.findAt(((-1.0, -1.0, 0.0), ))
mdb.models['Model-1'].DisplacementBC(name='PinY', createStepName='Load',
    region=corner, u2=0.0)

a.Surface(name='LoadSurface',
    side1Edges=a.instances['Plate-1'].edges.findAt(((1.0, 0.0, 0.0), )))
mdb.models['Model-1'].Pressure(name='Tension', createStepName='Load',
    region=a.surfaces['LoadSurface'], magnitude=1000000.0)

mdb.save()
