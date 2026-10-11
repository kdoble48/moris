"""Local Sphinx extension for the MORIS documentation site.

- Folds long Breathe "Subclassed by ..." lists into a collapsed <details> block:
  base classes such as fem::IWG have dozens of subclasses, printed as one
  paragraph ahead of the members.
- Registers placeholder Breathe directives when the Doxygen XML is absent
  (``make -C doc html-nodox``), so API pages render a note instead of failing.
"""

from docutils import nodes
from sphinx.util.docutils import SphinxDirective

SUBCLASS_FOLD_THRESHOLD = 8

BREATHE_DIRECTIVES = (
    "doxygenclass", "doxygenstruct", "doxygenfile", "doxygennamespace",
    "doxygenfunction", "doxygenenum", "doxygentypedef", "doxygendefine",
)


class folded(nodes.General, nodes.Element):
    """A block rendered as an HTML <details> element with a summary line."""


def visit_folded_html(self, node):
    self.body.append(f'<details class="moris-subclasses"><summary>{node["summary"]}</summary>')


def depart_folded_html(self, node):
    self.body.append("</details>")


def fold_subclass_lists(app, doctree):
    for paragraph in list(doctree.findall(nodes.paragraph)):
        text = paragraph.astext()
        if not text.startswith("Subclassed by "):
            continue
        count = text.count(",") + 1
        if count <= SUBCLASS_FOLD_THRESHOLD:
            continue
        block = folded(summary=f"Subclassed by {count} classes")
        paragraph.replace_self(block)
        block += paragraph


class _AnyOptions(dict):
    """Accept every option the real Breathe directives accept."""

    def __missing__(self, key):
        return lambda value: value


class DoxygenPlaceholder(SphinxDirective):
    required_arguments = 1
    final_argument_whitespace = True
    has_content = True
    option_spec = _AnyOptions()

    def run(self):
        text = f"{self.name}: {self.arguments[0]} (Doxygen XML not built; run `make -C doc html`)"
        return [nodes.note("", nodes.paragraph(text=text))]


def setup(app):
    app.add_config_value("moris_doxygen_available", True, "env")
    app.add_node(folded, html=(visit_folded_html, depart_folded_html))
    app.connect("doctree-read", fold_subclass_lists)
    app.connect("config-inited", _register_placeholders)
    return {"parallel_read_safe": True, "parallel_write_safe": True}


def _register_placeholders(app, config):
    if not config.moris_doxygen_available:
        for name in BREATHE_DIRECTIVES:
            app.add_directive(name, DoxygenPlaceholder, override=True)
