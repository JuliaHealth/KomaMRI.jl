"""
For using the Limit Struct of KomaMRICore when writing with MRIFiles
"""
insertNode(xs, paramName::String, param::KomaMRICore.Limit) = begin
    insertNode(xs, paramName, MRIFiles.Limit(param.minimum, param.maximum, param.center))
end
