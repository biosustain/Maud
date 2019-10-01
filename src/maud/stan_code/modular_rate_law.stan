real Tr_{{enz_id}} = p[{{enz}}]*p[{{Kcat1}}]
                *  {%- for su in substrate_list %}
                (m[{{su[0]}}]/p[{{su[1]}}])^(-1*{{su[2]}})
                {%- endfor %}
                - p[{{enz}}]*p[{{Kcat1}}]/p[{{Keq}}]
                *  {%- for su in substrate_list %}
                (p[{{su[1]}}])^{{su[2]}}
                {%- endfor %}
                *  {%- for pr in product_list %}
                (p[{{pr[1]}}])^{{pr[2]}}
                {%- endfor %}
                *  {%- for pr in product_list %}
                (m[{{pr[0]}}]/p[{{pr[1]}})^{{pr[2]}}
                {%- endfor %};

real Dr_{{enz_id}} = {%- for su in substrate_list %}
                    (1 + m[{{su[0]}}]/[p[{{su[1]}}])^(-1*{{su[2]}})
                {%- endfor %}
                + {%- for pr in product_list %}
                (1 + m[{{pr[0]}}]/p[{{pr[1]}})^{{pr[2]}}
                {%- endfor %}
                - 1;