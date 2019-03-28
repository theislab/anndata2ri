:github_url: {{ fullname | github_url }}

{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}
.. currentmodule:: {{ fullname }}

{% if classes -%}
Classes
-------

{%- for class in classes -%}
{%- if not class.startswith('_') %}
.. autoclass:: {{ class }}
   :members:
{%- endif -%}
{%- endfor -%}
{%- endif %}

{% if functions -%}
Functions
---------

{%- for function in functions -%}
{%- if not function.startswith('_') %}
.. autofunction:: {{ function }}
{%- endif -%}
{%- endfor -%}
{%- endif %}
