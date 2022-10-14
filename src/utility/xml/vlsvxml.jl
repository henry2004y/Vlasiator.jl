"""
    VLSVXML

Module of reading an XML-based VLSV file.

```julia-repl
julia> doc = readxml("test/init.vlsv");
julia> r = root(doc);
julia> ns = nodes(r);
```
"""
module VLSVXML

export
   # token
   Token, token_type,
   NUMBER, STRING, IDENT,
   UNKNOWN, ERROR, EOF,
   BEGIN_TAG, END_TAG, CLOSE_TAG, END_AND_CLOSE_TAG,
   # lexer
   Lexer, Token,
   lex, lex_end, scan_number, scan_string,
   ignore, ignore_whitespace,
   next_char, backup_char, peek_char, current_char,
   accept_char, accept_char_run,
   emit_token, lexeme,
   next_token,
   # parser
   Parser,
   peek_token, backup_token,
   peek_token_type,
   expect,
   # dom
   Node, ElementNode, OffsetNode, AttributeNode,
   nodename, iselement, istext, hasroot, hasnode,
   nodecontent,
   countnodes, countattributes,
   nodes, elements, attributes, eachattribute,
   root, setroot!,
   addchild!, addchildren!, addelement!,
   # XPath-like query API
   findall, findfirst,
   # xml_lexer
   lex_xml,
   # xml_parser
   parsexml, readxml,
   xmlparser, parse_node, parse_element, # Debug, remove later
   # show
   xml

include("token.jl")
include("lexer.jl")
include("parser.jl")
include("dom.jl")
include("xml_lexer.jl")
include("xml_parser.jl")
include("show.jl")

end