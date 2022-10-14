import Base: show

####### Show s-expression ###########
function show(io::IO, p::Node, depth::Integer=0)
   println(io)
   print(io, "  "^depth)
   print(io, "(", typeof(p), " \"")
   print(io, istext(p) ? nodecontent(p) : nodename(p))
   print(io, "\"")
   attrs = attributes(p)
   if !isempty(attrs)
      for a in attrs
         show(io, a, depth + 1) 
      end 
   end
    
   children = nodes(p)
   if isempty(children)
      print(io, ")")
   else
      for n in children
         show(io, n, depth + 1)
      end
      print(io, ")")
   end
end

function show(io::IO, n::AttributeNode, depth::Integer=0)
   println(io)
   print(io, "  "^depth)
   print(io, "(", typeof(n), " ")
   print(io, "\"", n.name, "\" = ", "\"", n.value, "\")")    
end