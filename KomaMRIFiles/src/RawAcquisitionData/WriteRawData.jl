"""
For using the Limit Struct of KomaMRIBase when writing with MRIFiles
"""
insertNode(xs, paramName::String, param::KomaMRIBase.Limit) = begin
    insertNode(xs, paramName, MRIFiles.Limit(param.minimum, param.maximum, param.center))
end
