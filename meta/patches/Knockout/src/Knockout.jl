module Knockout

using WebIO, Observables, JSExpr, JSON
import Observables: off, observe, AbstractObservable, ObservablePair

export knockout

#const knockout_js = joinpath(@__DIR__, "..", "assets", "knockout.js")
#const knockout_punches_js = joinpath(@__DIR__, "..", "assets", "knockout_punches.js")

"""
`knockout(template, data=Dict(), extra_js = js""; computed = [], methods = [])`

Create a Knockout scope, with HTML structure provided by `template` and filled with `data`.
# Arguments
- `template` the `Node` that acts as the template. See [Knockout syntax](http://knockoutjs.com/documentation/value-binding.html)
- `data` is either a dictionary or an array of `propertyName => value` pairs.
If a property's value is an observable, this function automatically sets up Julia -> JS communication.
To set up JS to Julia communication set up an event handler on `scope[propertyName]` (by calling `on(f, scope[propertyName])`)
_before_ rendering the scope.
You can specify that you want some knockout observable to be computed as a function of other observables,
e.g `knockout(...; computed = Dict(:fullName => @js function(){this.firstName() + ' ' + this.lastName()}))`.
You can pass functions that you want available in the Knockout scope as keyword arguments to
`knockout` E.g. `knockout(...; methods=Dict(:sayhello=>@js function(){ alert("hello!") }))`
"""
function knockout(template, data=Dict(), extra_js = js""; computed = [], methods = [])
    knockout_js = joinpath(Base.pkgdir(@__MODULE__), "assets", "knockout.js")
    knockout_punches_js = joinpath(Base.pkgdir(@__MODULE__), "assets", "knockout_punches.js")

    widget = Scope(imports=[
        "knockout" => knockout_js,
        "knockout_punches" => knockout_punches_js,
    ])
    widget(template)
    ko_data = Dict()
    watches = Dict()
    for (k, v) in data
        skey = string(k)
        (v isa ObservablePair) && (v = v.second)
        if isa(v, AbstractObservable)
            # associate the observable with the widget
            setobservable!(widget, skey, v)

            # forward updates from Julia to Knockoutjs
            onjs(v, @js function (val)
                if val != this.model[$skey]()
                    this.valueFromJulia[$skey] = true
                    this.model[$skey](val)
                end
            end)

            # forward updates from Knockoutjs to Julia
            watches[skey] = @js this[$skey].subscribe( function(val)
                if !this.valueFromJulia[$skey]
                    $v[] = val
                end
                this.valueFromJulia[$skey] = false
            end, self)
            ko_data[skey] = @js $v[]
        else
            ko_data[skey] = v
        end
    end

    methods_dict = Dict()
    for (k, f) in methods
        skey = string(k)
        methods_dict[skey] = @js this[$skey] = $f
    end

    computed_dict = Dict()
    for (k, f) in computed
        skey = string(k)
        computed_dict[skey] = @js this[$skey] = ko.computed($f, this)
    end

    on_import = js"""
    function (ko, koPunches) {
        ko.punches.enableAll();
        ko.bindingHandlers.numericValue = {
            init: function(element, valueAccessor, allBindings, data, context) {
                var stringified = ko.observable(ko.unwrap(valueAccessor()));
                stringified.subscribe(function(value) {
                    var val = parseFloat(value);
                    if (!isNaN(val)) {
                        valueAccessor()(val);
                    }
                });
                valueAccessor().subscribe(function(value) {
                    var str = JSON.stringify(value);
                    if ((str == "0") && (["-0", "-0."].indexOf(stringified()) >= 0))
                         return;
                     if (["null", ""].indexOf(str) >= 0)
                         return;
                    stringified(str);
                });
                ko.applyBindingsToNode(
                    element,
                    {
                        value: stringified,
                        valueUpdate: allBindings.get('valueUpdate'),
                    },
                    context,
                );
            }
        };
        var json_data = $ko_data;
        var self = this;
        function AppViewModel() {
            for (var key in json_data) {
                var el = json_data[key];
                this[key] = Array.isArray(el) ? ko.observableArray(el) : ko.observable(el);
            }
            $(dict2js(methods_dict))
            $(dict2js(computed_dict))
            $(dict2js(watches))
            $extra_js
        }
        self.model = new AppViewModel();
        self.valueFromJulia = {};
        for (var key in json_data) {
            self.valueFromJulia[key] = false;
        }
        ko.applyBindings(self.model, self.dom);
    }
    """
    onimport(widget, on_import)
    widget
end

function dict2js(d::AbstractDict)
    isempty(d) ? js"" : js"$(values(d)...,)"
end

isnumeric(x) = false
isnumeric(x::Number) = true
isnumeric(x::Bool) = false
isnumeric(x::AbstractObservable) = isnumeric(x[])

js_lambda(s::String) = "function (){$s}"

end # module
