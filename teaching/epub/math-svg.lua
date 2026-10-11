-- Render every TeX formula to an SVG image with xelatex + dvisvgm.
-- xelatex (not latex) because some formulas contain Korean text.
-- All formulas are typeset in one latex run (one preview page each), so the
-- build takes seconds. Each image is sized in em from its TeX box, and its
-- depth sets vertical-align, so inline math sits on the text baseline.

local preamble = [[
\documentclass[12pt]{article}
\usepackage{amsmath,amssymb}
\usepackage{kotex}
\setmainhangulfont{NanumMyeongjo}
\usepackage[active,tightpage,displaymath,textmath]{preview}
\setlength\PreviewBorder{0pt}
\newcommand{\logit}{\operatorname{logit}}
\newcommand{\Var}{\operatorname{Var}}
\begin{document}
]]

local function collect(doc)
  local items, index = {}, {}
  doc:walk({
    Math = function(m)
      local key = m.mathtype .. "\0" .. m.text
      if not index[key] then
        items[#items + 1] = m
        index[key] = #items
      end
    end,
  })
  return items, index
end

function Pandoc(doc)
  local items, index = collect(doc)
  if #items == 0 then return doc end
  local dir = os.getenv("MATH_DIR") or "math"
  os.execute("mkdir -p " .. dir)
  local body = {}
  for _, m in ipairs(items) do
    if m.mathtype == "DisplayMath" then
      body[#body + 1] = "\\[" .. m.text .. "\\]\n"
    else
      -- Measure the depth outside the preview so the image can sit on the
      -- baseline; the preview page itself is exactly the formula's box.
      body[#body + 1] = "\\begin{nopreview}\\sbox0{$" .. m.text .. "$}\\typeout{EQDEPTH=\\the\\dp0}\\end{nopreview}\n$" .. m.text .. "$\n\n"
    end
  end
  local tex = io.open(dir .. "/all.tex", "w")
  tex:write(preamble, table.concat(body, "\n"), "\\end{document}\n")
  tex:close()
  assert(os.execute("cd " .. dir .. " && xelatex -no-pdf -interaction=nonstopmode -halt-on-error all.tex > latex.log 2>&1"),
         "latex failed; see " .. dir .. "/latex.log")
  local handle = io.popen("cd " .. dir .. " && dvisvgm --page=1- --bbox=papersize --no-fonts=1 -o 'eq-%p.svg' all.xdv 2>&1")
  local log = handle:read("a")
  handle:close()
  -- dvisvgm prints "graphic size: <w>pt x <h>pt" per page; the TeX log holds
  -- the inline depths, in formula order.
  local depths = {}
  for d in io.open(dir .. "/all.log"):read("a"):gmatch("EQDEPTH=([%d%.]+)pt") do
    depths[#depths + 1] = tonumber(d)
  end
  local sizes, inline = {}, 0
  for w, h in log:gmatch("graphic size: ([%d%.]+)pt x ([%d%.]+)pt") do
    local d = 0
    if items[#sizes + 1] and items[#sizes + 1].mathtype == "InlineMath" then
      inline = inline + 1
      d = depths[inline] or 0
    end
    sizes[#sizes + 1] = { tonumber(w), tonumber(h), d }
  end
  assert(#sizes == #items, ("expected %d formulas, dvisvgm reported %d"):format(#items, #sizes))
  local digits = #tostring(#items)
  return doc:walk({
    Math = function(m)
      local i = index[m.mathtype .. "\0" .. m.text]
      local w, h, d = table.unpack(sizes[i])
      local file = ("%s/eq-%0" .. digits .. "d.svg"):format(dir, i)
      local em = 12  -- formulas are typeset at 12pt
      local style = ("height:%.3fem;vertical-align:%.3fem"):format(h / em, -d / em)
      local img = pandoc.Image({ pandoc.Str(m.text) }, file, "", { class = "math " .. m.mathtype, style = style })
      if m.mathtype == "DisplayMath" then
        return pandoc.Span({ img }, { class = "display-math" })
      end
      return img
    end,
  })
end
