# Cantilever beam under tip load — classic Abaqus tutorial script.
# SI units (m, Pa, N).
from abaqus import *
from abaqusConstants import *
import mesh

s = mdb.models['Model-1'].ConstrainedSketch(name='profile', sheetSize=2.0)
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0.0, 0.0), point2=(1.0, 0.1))
p = mdb.models['Model-1'].Part(name='Beam',
    dimensionality=TWO_D_PLANAR,
    type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)
s.unsetPrimaryObject()

mdb.models['Model-1'].Material(name='Steel')
mdb.models['Model-1'].materials['Steel'].Elastic(table=((210000000000.0, 0.3), ))
mdb.models['Model-1'].materials['Steel'].Density(table=((7850.0,), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='SteelSection',
    material='Steel', thickness=0.01)
p.SectionAssignment(region=p.Set(faces=p.faces, name='AllFaces'),
    sectionName='SteelSection')

p.seedPart(size=0.025, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

a = mdb.models['Model-1'].rootAssembly
a.Instance(name='Beam-1', part=p, dependent=ON)

mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial',
    description='tip load')

fixRegion = a.instances['Beam-1'].edges.findAt(((0.0, 0.05, 0.0), ))
mdb.models['Model-1'].DisplacementBC(name='Fix', createStepName='Step-1',
    region=fixRegion, u1=0.0, u2=0.0)

tipRegion = a.instances['Beam-1'].vertices.findAt(((1.0, 0.05, 0.0), ))
mdb.models['Model-1'].ConcentratedForce(name='TipLoad', createStepName='Step-1',
    region=tipRegion, cf2=-1000.0)

mdb.save()
