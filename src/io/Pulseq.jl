# # add_format(format"Pulseq", "# Pulseq sequence file", ".seq", [:Koma=>UUID("a9056882-1c2f-47b4-848d-dcb4a04f1994")])
# module Pulseq

# using FileIO

# # Again, this is a *private* `load` function, do not extend `FileIO.load`!
# function load(f::File{format"Pulseq"})
#     open(f) do s
#         skipmagic(s)  # skip over the magic bytes
#         # You can just call the `load(::Stream)` method below...
#         ret = load(s)
#         # ...or implement everything here instead
#     end
# end

# # You can support streams and add keywords:
# function load(s::Stream{format"Pulseq"}; keywords...)
#     # s is already positioned after the magic bytes
#     # Do the stuff to read a PNG file
#     chunklength = read(s, UInt32)
#     ...
# end

# function save(f::File{format"Pulseq"}, data)
#     open(f, "w") do s
#         # Don't forget to write the magic bytes!
#         write(s, magic(format"PNG"))
#         # Do the rest of the stuff needed to save in PNG format
#     end
# end

# end 