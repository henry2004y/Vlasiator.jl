tagstrip(tag::AbstractString) = strip(tag, ['<', '>', '/'])

function match_closing_tag(n::ElementNode, t::Token)
   close_tagname = tagstrip(t.lexeme)
   if n.name != close_tagname
      error("opening tag '$(n.name)' does not match closing tag '$close_tagname'")
   end
end

function add_parsed_nodes!(parser::Parser, parent::ElementNode)
   t = peek_token(parser)
   while t.kind == BEGIN_TAG || t.kind == TEXT
      addchild!(parent, parse_node(parser)) #!!!
      t = peek_token(parser)
   end
end

function add_parsed_attribute_nodes!(parser::Parser, parent::ElementNode)
   t = peek_token(parser)
   while t.kind == IDENT
      t = expect(parser, IDENT)
      name = t.lexeme
      next_token(parser)
      t = expect(parser, STRING)
      value = t.lexeme
      push!(parent.attributes, AttributeNode(name, value))
      t = peek_token(parser)
   end
end

function parse_text(parser::Parser)
   t = expect(parser, TEXT)
   parse(Int, t.lexeme) |> OffsetNode
end

# Most time-consuming function!
function parse_element(parser::Parser)
   t = peek_token(parser)
   expect(parser, BEGIN_TAG)
   tagname = tagstrip(t.lexeme)
   n = ElementNode(tagname)

   t = peek_token(parser)
   while t.kind in [END_TAG, IDENT]
      if t.kind == END_TAG
         expect(parser, END_TAG)
         add_parsed_nodes!(parser, n) #!!!
      elseif t.kind == IDENT
         add_parsed_attribute_nodes!(parser, n) #!!!
      end
      t = peek_token(parser)
   end

   t = next_token(parser)
   if t.kind == CLOSE_TAG
      match_closing_tag(n, t)
   end

   if t.kind ∉ [CLOSE_TAG, END_AND_CLOSE_TAG]
      error("Element node needs to end with /> or </$(n.name)> not '$t'")
   end

   return n
end

function parse_node(parser::Parser)
   t = peek_token(parser)
   if t.kind == BEGIN_TAG
      parse_element(parser)::Node #!!!
   elseif t.kind == TEXT
      parse_text(parser)::Node
   else
      error("Had not expected token '$t' while looking for start of new XML node")
   end
end

xmlparser(s::AbstractString) = Parser(lex_xml(s))

"""
    parsexml(xmlstring::AbstractString; ignore_declaration=false)
    
Parse a text string containing XML, and return a vector of ElementNode.
"""
function parsexml(s::AbstractString)
   l = lex_xml(s)
   p = Parser(l)
   parse_element(p)
end

"""
    readxml(stream::IO=stdin) -> Vector{ElementNode}
    readxml(filename::AbstractString) -> Vector{ElementNode}
    
Read and parse XML from an I/O stream or file.
"""
function readxml(stream::IO=stdin)
   text = read(stream, String)
   parsexml(text)
end

readxml(filename::AbstractString) = open(readxml, filename)