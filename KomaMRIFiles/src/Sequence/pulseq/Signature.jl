function supported_signature_digest(algorithm::AbstractString, payload::Vector{UInt8})
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

function parse_signature_section(section::AbstractString)
    io = IOBuffer(section)
    readline(io) # consume "[SIGNATURE]" header
    signature_type = nothing
    signature_hash = nothing
    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue
        startswith(line, '#') && continue
        parts = split(line)
        isempty(parts) && continue
        key = parts[1]
        value = join(parts[2:end], " ")
        if key == "Type"
            signature_type = strip(value)
        elseif key == "Hash"
            signature_hash = lowercase(replace(strip(value), " " => ""))
        end
    end
    return isnothing(signature_type) && isnothing(signature_hash) ? nothing : (type = signature_type, hash = signature_hash)
end

function verify_signature!(filename::String, signature::Nothing; pulseq_version::VersionNumber=v"1.4.0")
    @warn "Pulseq [SIGNATURE] section is missing; skipping verification."
    return nothing
end

function verify_signature!(filename::String, signature::NamedTuple; pulseq_version::VersionNumber=v"1.4.0")
    # Don't error, just warn - the file can still be used without the signature or with it being incorrect
    sig_type = signature.type
    sig_hash = signature.hash
    # Read file as bytes to avoid any line ending normalization
    file_bytes = read(filename)
    # Find [SIGNATURE] in bytes
    sig_marker = b"[SIGNATURE]"
    sig_pos = findfirst(sig_marker, file_bytes)
    isnothing(sig_pos) && begin
        @warn "Signature section expected but not found when verifying Pulseq file." filename
        return
    end
    sig_start = first(sig_pos)
    # Extract payload based on version
    # For Pulseq < 1.4.0 (e.g., JEMRIS), the newline before [SIGNATURE] is part of the payload
    # For Pulseq >= 1.4.0, the newline before [SIGNATURE] is part of the signature and should be excluded
    # However, different implementations may handle this differently, so we try multiple approaches
    payload_end = sig_start - 1
    expected_hash = lowercase(replace(sig_hash, " " => ""))
    if pulseq_version < v"1.4.0"
        # Include the newline before [SIGNATURE] in the payload (JEMRIS format)
        payload_bytes = payload_end > 0 ? file_bytes[1:payload_end] : UInt8[]
        computed_hash = supported_signature_digest(sig_type, payload_bytes)
    else
        # For version >= 1.4.0, try multiple approaches to handle different implementations
        # 1. Exclude only the last newline (spec-compliant)
        if payload_end > 0 && file_bytes[payload_end] in (UInt8('\n'), UInt8('\r'))
            payload_bytes_excluding = file_bytes[1:(payload_end - 1)]
        else
            payload_bytes_excluding = payload_end > 0 ? file_bytes[1:payload_end] : UInt8[]
        end
        # 2. Include the last newline (some implementations like MATLAB)
        payload_bytes_including = file_bytes[1:payload_end]
        # 3. Exclude all consecutive newlines before [SIGNATURE] (some implementations)
        last_non_nl = payload_end
        while last_non_nl > 0 && file_bytes[last_non_nl] in (UInt8('\n'), UInt8('\r'))
            last_non_nl -= 1
        end
        payload_bytes_no_nl = last_non_nl > 0 ? file_bytes[1:last_non_nl] : UInt8[]
        # Try all approaches and use the one that matches
        computed_hash_excluding = supported_signature_digest(sig_type, payload_bytes_excluding)
        computed_hash_including = supported_signature_digest(sig_type, payload_bytes_including)
        computed_hash_no_nl = supported_signature_digest(sig_type, payload_bytes_no_nl)
        if computed_hash_excluding == expected_hash
            computed_hash = computed_hash_excluding
            payload_bytes = payload_bytes_excluding
        elseif computed_hash_including == expected_hash
            computed_hash = computed_hash_including
            payload_bytes = payload_bytes_including
        elseif computed_hash_no_nl == expected_hash
            computed_hash = computed_hash_no_nl
            payload_bytes = payload_bytes_no_nl
        else
            # None matches, use the spec-compliant one for error message
            computed_hash = computed_hash_excluding
            payload_bytes = payload_bytes_excluding
        end
    end
    # Parse signature section for comparison (read as string for parsing)
    file_text = String(file_bytes)
    sig_section_start = sig_start
    sig_section = file_text[sig_section_start:end]
    parsed = parse_signature_section(sig_section)
    isnothing(parsed) && begin
        @warn "Failed to parse Pulseq signature section; skipping verification." filename
        return
    end
    parsed_type = parsed.type
    parsed_hash = parsed.hash
    if !isnothing(parsed_type) && !isempty(parsed_type) && lowercase(parsed_type) != lowercase(sig_type)
        @warn "Pulseq signature Type mismatch between parser and metadata. Using '$sig_type'." filename parsed_type
    end
    if !isnothing(parsed_hash) && !isempty(parsed_hash)
        parsed_hash_norm = lowercase(replace(parsed_hash, " " => ""))
        if parsed_hash_norm != lowercase(replace(sig_hash, " " => ""))
            @warn "Pulseq signature Hash mismatch between parser and metadata. Using metadata value." filename
        end
    end
    if computed_hash != expected_hash
        @warn "Pulseq signature verification failed for $(basename(filename)). Expected $(expected_hash), computed $(computed_hash). The file may have been modified or generated with a different implementation." filename
    end
end