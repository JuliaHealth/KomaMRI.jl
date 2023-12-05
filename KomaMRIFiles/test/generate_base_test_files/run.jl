# This script saves the following files:
# * spiral.jld2
# * epi.jld2
# * radial_JEMRIS.jld2
# * epi_100x100_TE100_FOV230.jld2
# * sphere_chemical_shift.jld2
# * jemris_signals_epi_sphere_cs.jld2
# * Koma_signal.jld2
# in the directory KomaMRIBase/test/test_files/
#
# Note that the saved phantom sphere_chemical_shift.jld2 doesn't have stored the information
# of the motion functions

using JLD2, HDF5, KomaMRIBase, KomaMRIFiles
path = @__DIR__
path_mrifiles = path*"/../test_files"
path_mribase = path*"/../../../KomaMRIBase/test/test_files"
mkpath(path_mribase)

# Save sequence structs to files
seq = read_seq(path_mrifiles*"/spiral.seq")
save_object(path_mribase*"/spiral.jld2", seq)
seq = read_seq(path_mrifiles*"/epi.seq")
save_object(path_mribase*"/epi.jld2", seq)
seq = read_seq(path_mrifiles*"/radial_JEMRIS.seq")
save_object(path_mribase*"/radial_JEMRIS.jld2", seq)
seq = read_seq(path_mrifiles*"/epi_100x100_TE100_FOV230.seq")
save_object(path_mribase*"/epi_100x100_TE100_FOV230.jld2", seq)

# Save phantom structs to files
# (there is a problem when reading the functions, so we don't save them)
function phantom_to_namedtuple(obj::Phantom)
    fields = fieldnames(typeof(obj))[1:end-3]
    values = [getfield(obj, field) for field in fields]
    return NamedTuple{fields}(values)
end
obj = read_phantom_jemris(path_mrifiles*"/sphere_chemical_shift.h5")
tuple_obj = phantom_to_namedtuple(obj)
save_object(path_mribase*"/sphere_chemical_shift.jld2", tuple_obj)

# Save jemris signal
sig_jemris = h5open(path_mrifiles*"/jemris_signals_epi_sphere_cs.h5")["/signal/channels/00"]
sig_jemris = sig_jemris[1,:] + 1im*sig_jemris[2,:]
sig_jemris = sig_jemris[:]
save_object(path_mribase*"/jemris_signals_epi_sphere_cs.jld2", sig_jemris)

# Save Koma_signal
raw = RawAcquisitionData(ISMRMRDFile(path_mrifiles*"/Koma_signal.mrd"))
save_object(path_mribase*"/Koma_signal.jld2", raw)
