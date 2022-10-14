# NOTE: The interface has been modeled as far as possible on the EzXML.jl XML
# parser. This is so EzXML can be a drop in replacement. The benefits of this
# parser over EzXML is that it has no dependencies. It is all pure Julia code.
# The downside is that it just supports the most common XML features
import Base: getindex, haskey, findfirst

"All nodes in XML DOM is some type of Node."
abstract type Node end

"""
    AttributeNode

Represents an attribute inside a tag.

# Example
Here `class` and `name` are examples of attributs belonging to parent node `widget`.

    <widget class="QCheckBox" name="checkBox">
"""
struct AttributeNode <: Node
   name::String
   value::String
end

"XML Node which can contain attributes and child nodes"
struct ElementNode <: Node
   name::String
   attributes::Vector{AttributeNode}
   children::Vector{Node}
end

function ElementNode(name::AbstractString)
   ElementNode(name, AttributeNode[], Node[])
end

"Offset for each parameter/variable."
struct OffsetNode <: Node
   content::Int
end

"Creates an element node named `name` with a text node containing text `value`"
function ElementNode(name::AbstractString, value::AbstractString)
   ElementNode(name, AttributeNode[], Node[OffsetNode(value)])
end

"Element node with children `nodes`"
function ElementNode(name::AbstractString, nodes::Array{T}) where T <: Node
   ElementNode(name, AttributeNode[], nodes)
end

"""
Element node with attributes given like
`ElementNode("widget", ["class"=>class, "name"=>name])`
"""
function ElementNode(name::AbstractString, attributes::Vector{Pair{String, String}})
   ElementNode(name, [AttributeNode(name, value) for (name, value) in attributes], Node[])
end

## Public API

function getindex(n::ElementNode, key::AbstractString)
   for m in n.attributes
      if m.name == key
         return m.value
      end
   end
   error("No attribute with key $key exist")
end

getindex(n::ElementNode, i::Integer) = n.children[i]


"Get all child nodes under node `n`"
nodes(n::ElementNode)::Vector{VLSVXML.ElementNode} = n.children

"Get all elements under node `n`"
elements(n::ElementNode) = filter(iselement, nodes(n))

"Get an array of attributes under node `n`"
attributes(n::ElementNode) = n.attributes

"Gets a dictionary of attributes meant to use in a for loop for iteration"
eachattribute(n::ElementNode) = n.attributes

"For an XML tag looking like `<foo>bar</foo>` the `nodename` would be foo"
nodename(n::OffsetNode) = "offset"
nodename(n::ElementNode) = n.name

"Check if node is an element node. Element nodes can contain child nodes"
iselement(n::Node) = false
iselement(n::ElementNode) = true

"Number of child nodes"
countnodes(n::Node) = 0
countnodes(n::ElementNode) = length(n.children)

"Number of attributes. Typically onlye element nodes have attributes"
countattributes(n::Node) = 0
countattributes(n::ElementNode) = length(n.attributes)

"Get content of all text nodes under `n`"
nodecontent(n::OffsetNode) = n.content

"""
    hasnode(node)
Return if `node` has a child node.
"""
hasnode(n::Node) = false
hasnode(n::ElementNode) = !isempty(n.children)

function addelement!(parent::Node, name::AbstractString)
   child = ElementNode(name)
   addchild!(parent, child)
   child
end

"Add `child` node to `parent` node"
addchild!(parent::ElementNode, child::Node) = push!(parent.children, child)

## XPath-like API
# The standard XPath starts with "//". Here we do not follow this rule but keep the name.

"""
    haskey(node, attribute) -> Bool
Check if XML node has a particular attribute. E.g.`haskey(n, "foobar")`
would return `true` for `<egg foobar="spam"/>` but `false` for `<foobar egg="spam"/>`
"""
haskey(n::Node, key::AbstractString) = any(m->m.name == key, n.attributes)

function Base.findall(xpath::AbstractString, node::ElementNode)
   [node.children[i] for i in eachindex(node.children) if node.children[i].name == xpath]
end

function Base.findfirst(xpath::AbstractString, node::ElementNode)
   ns = findall(xpath, node)
   isempty(ns) ? nothing : first(ns)
end