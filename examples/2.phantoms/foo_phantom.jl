using HDF5

# Create HDF5 phantom file
fid = h5open("foo.phantom","w")

# Root attributes
attributes(fid)["Version"] = "1.0"
attributes(fid)["Name"] = "Displacement"
attributes(fid)["Ns"] = 2
attributes(fid)["Dims"] = 2
attributes(fid)["Dynamic"] = 1     # 0=False, 1=True

# Spin initial positions
pos = create_group(fid,"position")
pos["x"] = [0, 0]
pos["y"] = [-0.05, 0.05]

# Contrast (Rho, T1, T2, Deltaw)
contrast = create_group(fid,"contrast")

rho = create_group(contrast,"rho")
attributes(rho)["type"] = "Explicit"
rho["values"] = [1,1]

T1 = create_group(contrast,"T1")
attributes(T1)["type"] = "Indexed"
T1["values"] = [1,2]
T1["table"] = [3,4]
attributes(T1)["N"] = length(T1["table"])

T2 = create_group(contrast,"T2")
attributes(T2)["type"] = "Explicit"
T2["values"] = [1,1]

Deltaw = create_group(contrast,"Deltaw")
attributes(Deltaw)["type"] = "Explicit"
Deltaw["values"] = [0 0]

# Motion
motion = create_group(fid,"motion")

segments = create_group(motion, "segments")
attributes(segments)["N"] = 1
attributes(segments)["K"] = 4
segments["dur"] = [0.1]


motionx = create_group(motion,"motionx")
attributes(motionx)["type"] = "Explicit"
motionx["values"] = [-0.1 0 0.1;
                     0.1 0 -0.1]

close(fid)
