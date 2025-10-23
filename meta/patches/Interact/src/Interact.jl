module Interact

using Reexport

@reexport using InteractBase
import InteractBase: notifications
import Widgets: Widget, @layout, @nodeps
import Observables: @on, @map!, @map

@reexport using OrderedCollections
@reexport using Observables
@reexport using Knockout
@reexport using CSSUtil
@reexport using WebIO
@reexport using Widgets

struct Bulma<:InteractBase.WidgetTheme; end

# const notebookdir = joinpath(@__DIR__, "..", "doc", "notebooks")

# const bulma_css = joinpath(@__DIR__, "..", "assets", "bulma.min.css")
# const bulma_confined_css = joinpath(@__DIR__, "..", "assets", "bulma_confined.min.css")

function InteractBase.libraries(::Bulma)
    bulma_css = joinpath(Base.pkgdir(@__MODULE__), "assets", "bulma.min.css")
    bulma_confined_css = joinpath(Base.pkgdir(@__MODULE__), "assets", "bulma_confined.min.css")
    
    bulmalib = InteractBase.isijulia() ? bulma_confined_css : bulma_css
    vcat(InteractBase.font_awesome, InteractBase.style_css, bulmalib)
end

function __init__()

    root = joinpath(dirname(Base.pkgdir(InteractBase)), "Interact") # an ugly hack, again!
    global notebookdir = joinpath(root, "doc", "notebooks")

    InteractBase.registertheme!(:bulma, Bulma())
    settheme!(Bulma())
    nothing
end

end
