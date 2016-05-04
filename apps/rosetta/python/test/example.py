from pysetta import devel
from pysetta.core.pose import Pose
from pysetta.core.import_pose import pose_from_file

devel.init("dummy -database /work/sheffler/ssd/rosetta/librosetta/database".split())
p = Pose()
p.dump_pdb("test.pdb")
p = pose_from_file("/work/sheffler/1ffw_native.pdb")
p.dump_pdb("test2.pdb")
