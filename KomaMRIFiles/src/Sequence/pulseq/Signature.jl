const PulseqSignature = NamedTuple{(:type, :hash),Tuple{String,String}}

function supported_signature_digest(algorithm, payload)
    alg = lowercase(strip(algorithm))
    digest = if alg == "md5"
        md5(payload)
    elseif alg == "sha1"
        sha1(payload)
    elseif alg in ("sha2", "sha256")
        sha256(payload)
    else
        throw(ArgumentError("Unsupported signature algorithm '$algorithm'. Supported algorithms: md5, sha1, sha256."))
    end
    return lowercase(bytes2hex(digest))
end

function verify_signature!(filename, signature::Nothing; pulseq_version=v"1.4.0")
    @warn "Pulseq [SIGNATURE] section is missing; skipping verification."
    return nothing
end

function signature_payload_candidates(file_bytes, sig_start, pulseq_version)
    payload_end = sig_start - 1
    payload = payload_end > 0 ? file_bytes[1:payload_end] : UInt8[]
    pulseq_version < v"1.4.0" && return (payload,)
    payload_without_last_newline =
        payload_end > 0 && file_bytes[payload_end] in (UInt8('\n'), UInt8('\r')) ?
        file_bytes[1:(payload_end - 1)] :
        payload
    last_non_newline = payload_end
    while last_non_newline > 0 && file_bytes[last_non_newline] in (UInt8('\n'), UInt8('\r'))
        last_non_newline -= 1
    end
    payload_without_trailing_newlines = last_non_newline > 0 ? file_bytes[1:last_non_newline] : UInt8[]
    return (payload_without_last_newline, payload, payload_without_trailing_newlines)
end

function verify_signature!(filename, signature::PulseqSignature; pulseq_version=v"1.4.0")
    file_bytes = read(filename)
    sig_marker = b"[SIGNATURE]"
    sig_pos = findfirst(sig_marker, file_bytes)
    isnothing(sig_pos) && begin
        @warn "Signature section expected but not found when verifying Pulseq file." filename
        return
    end
    expected_hash = lowercase(replace(signature.hash, " " => ""))
    computed_hash = ""
    for payload in signature_payload_candidates(file_bytes, first(sig_pos), pulseq_version)
        computed_hash = supported_signature_digest(signature.type, payload)
        computed_hash == expected_hash && return nothing
    end
    if computed_hash != expected_hash
        @warn "Pulseq signature verification failed for $(basename(filename)). Expected $(expected_hash), computed $(computed_hash). The file may have been modified or generated with a different implementation." filename
    end
    return nothing
end
