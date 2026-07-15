setup_scanner() = Scanner()

setup_sequence(sys::Scanner) = PulseDesigner.EPI_example(; sys)

function setup_phantom(; phantom_mode="2D")
    obj = phantom_mode == "3D" ? brain_phantom3D() : brain_phantom2D()
    obj.Δw .= 0
    return obj
end

function setup_raw()
    return RawAcquisitionData(
        Dict(
            "systemVendor" => "",
            "encodedSize" => [2, 2, 1],
            "reconSize" => [2, 2, 1],
            "number_of_samples" => 4,
            "encodedFOV" => [100.0, 100.0, 1],
            "trajectory" => "other",
        ),
        [
            KomaMRICore.Profile(
                AcquisitionHeader(; trajectory_dimensions=2, sample_time_us=1),
                [0.0 0.0 1 1; 0 1 1 1] ./ 2,
                reshape([0.0; 0im; 0; 0], 4, 1),
            ),
        ],
    )
end

const sys_ui = Observable{Scanner}(Scanner())
const seq_ui = Observable{Sequence}(Sequence())
const obj_ui = Observable{Phantom}(Phantom(x=[0.0]))
const physio_ui = Observable{AbstractPhysioSignal}(NoPhysioSignal())
const raw_ui = Observable{RawAcquisitionData}(setup_raw())
const img_ui = Observable{AbstractArray{<:Complex}}([0.0im 0.; 0. 0.])

default_physio_signal(seq) = has_trigger(seq) ? CardiacSignal(; heart_rate=1) : NoPhysioSignal()

const MAX_UI_FILENAME_CHARS = 32
const MAX_TOAST_FILENAME_CHARS = 20
display_filename(name, limit=MAX_UI_FILENAME_CHARS) =
    length(name) > limit ? string(first(name, limit - 3), "...") : name
