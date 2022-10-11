import Base: ==

@enum TokenType begin
   NUMBER
   STRING
   IDENT             # Generic
   UNKNOWN
   ERROR
   EOF               # Control
   # XML
   BEGIN_TAG         # <tag
   END_TAG           # >
   CLOSE_TAG         # </tag>
   END_AND_CLOSE_TAG # />
   TEXT              # Text such as "bar" inside tags: <foo>bar</foo>
   EQUAL = Int('=')
end

"""
A string of code is turned into an array of Tokens by the lexer. Each symbol or word in
code is represented as a Token.
"""
struct Token
   kind::TokenType
   lexeme::String
   function Token(kind::TokenType, s::AbstractString)
      # Don't want to store the quotes on quoted strings.
      if length(s) > 1
         new(kind, strip(s, '"'))
      else
         new(kind, s)
      end
   end
end

function Token(kind::TokenType)
   ch = Char(kind)
   if ch in "{}(),=;"
      Token(kind, string(ch))
   else
      Token(kind, "")
   end
end

==(t1::Token, t2::Token) = t1.kind == t2.kind && t1.lexeme == t2.lexeme

token_type(token::Token) = token.kind