"""
For using the Limit Struct of KomaMRICore when writing with MRIFiles
"""
function MRIFiles.insertNode(xs, paramName::String, param::KomaMRICore.Limit)
    xsp = MRIFiles.new_child(xs, paramName)
    MRIFiles.insertNode(xsp, "minimum", param.minimum)
    MRIFiles.insertNode(xsp, "maximum", param.maximum)
    MRIFiles.insertNode(xsp, "center", param.center)
end
