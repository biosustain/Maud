real Tr_{{enz_id}} = p[{{enz}}]*p[{{Kcat1}}]
                *  {%- for su in substrate_list %} ({{su[0]}}/p[{{su[1]}}])^(-1*{{su[2]}}) {{"*" if not loop.last}}
                {%- endfor %}
                - p[{{enz}}]*p[{{Kcat1}}]/p[{{Keq}}]
                *  {%- for su in substrate_list %} (p[{{su[1]}}])^{{su[2]}} {{"*" if not loop.last}}
                {%- endfor %}
                *  {%- for pr in product_list %} (p[{{pr[1]}}])^{{pr[2]}} {{"*" if not loop.last}}
                {%- endfor %}
                *  {%- for pr in product_list %} ({{pr[0]}}/p[{{pr[1]}}])^{{pr[2]}} {{"*" if not loop.last}}
                {%- endfor %};

real Dr_{{enz_id}} = {%- for su in substrate_list %} (1 + {{su[0]}}/p[{{su[1]}}])^(-1*{{su[2]}}) {{"*" if not loop.last}}
                {%- endfor %}
                + {%- for pr in product_list %} (1 + {{pr[0]}}/p[{{pr[1]}}])^{{pr[2]}} {{"*" if not loop.last}}
                {%- endfor %}
                - 1;

{% if competitive_inhibitor_list %}
real Dr_reg_{{enz_id}} = {%- for ci in competitive_inhibitor_list %} ({{ci[0]}}/p[{{ci[1]}}]) {{"*" if not loop.last}}
    {%- endfor %};

{% else %}
real Dr_reg_{{enz_id}} = 0;
    {%- endif %}