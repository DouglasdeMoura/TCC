/*
 * Author: Douglas Araujo de Moura (douglas.moura [at] constrinew.com.br)
 * Author URI: https://engenharialivre.com/
 * License: GNU General Public License v2 or later
 * License URI: http://www.gnu.org/licenses/gpl-2.0.html
 * Repository: https://github.com/DouglasdeMoura/TCC
 * Last modified: 2017-11-14
*/

function x = resultados(b, h, bf, hf, g, q, l, fck)
    compressao_maxima = 0.6 * fck * 1000

    [Ac, I, Wsup, Winf] = caractGeometricas(b, h, bf, hf)
    [Mgmax, Mqmax, sigma_cg_inf, sigma_cg_sup, sigma_cq_inf, sigma_cq_sup, Npd] = esforcosSolicitantes(g, q, l, Wsup, Winf)

    Np = ( - sigma_cg_inf - sigma_cq_inf ) * Ac

    sigma_c_min = Npd + sigma_cg_sup + sigma_cq_sup

    if abs(sigma_c_min) < abs(compressao_maxima) then
        v1 = 1
    else
        v1 = 0
        //warning("A tensão na fibra superior (|" + string(sigma_c_min/1000)  + "| MPa) é maior que a máxima permitida (|" + string(compressao_maxima/1000)  + "| MPa).")
    end

    if abs(Npd) < abs(compressao_maxima) then
        v2 = 1
    else
        v2 = 0
        //warning("A compressão na seção dos apoios (|" + string(Npd/1000)  + "| MPa) é maior que a máxima permitida (|" + string(compressao_maxima/1000)  + "| MPa).")
    end

    // 1. Normal de protensão (MPa)
    // 2. fck do concreto (MPa)
    // 3. Compressão máxima permitida (MPa)
    // 4. Compressão mínima no concreto (MPa)
    // 5. Taxa de aproveitamento da compressão no concreto (%)
    // 6. Vão teórico (m)
    // 7. Status da primeira validação (boolean)
    // 8. Status da segunda validação (boolean)
    // 9. Área da seção transversal
    // 10. Inércia da seção
    // 11. W_sup
    // 12. W_inf
    x = [(Npd/1000) fck (abs(compressao_maxima)/1000) ((abs(sigma_c_min))/1000) ((abs(sigma_c_min) / abs(compressao_maxima)) * 100) l v1 v2 Ac I Wsup Winf]

endfunction

function [Ac, I, Wsup, Winf] = caractGeometricas(b, h, bf, hf)
    if bf == 0 & hf == 0 then
        // Seção retangular
        Ac = b * h
        I  = (b * (h^3)) / 12
        y = h / 2
        Wsup = I / y
        Winf = I / y
    else
        // Seção T
        h     = (h - hf)
        A1    = bf * hf
        A2    = b  * h
        Ac    = A1 + A2
        Msx   = (A1 * (h + (hf / 2))) + (A2 * (h / 2))
        yg    = Msx / Ac
        I_xg1 = ((bf * (hf^3)) / 12) + (bf * hf * (((h + (hf / 2)) - yg)^2))
        I_xg2 = ((b * (h^3)) / 12) + (b * h * ((yg - (h / 2))^2))
        I     = I_xg1 + I_xg2
        ysup  = (h + hf) - yg
        yinf  = yg
        Wsup  = I / ysup
        Winf  = I / yinf
    end
endfunction

function [Mgmax, Mqmax, sigma_cg_inf, sigma_cg_sup, sigma_cq_inf, sigma_cq_sup, Npd] = esforcosSolicitantes(g, q, l, Wsup, Winf)
    j = (l^2 / 8)
    Mgmax = g * j
    Mqmax = q * j
    
    sigma_cg_inf = + Mgmax / Winf
    sigma_cg_sup = - Mgmax / Wsup
    sigma_cq_inf = + Mqmax / Winf
    sigma_cq_sup = - Mqmax / Wsup
    Npd = ( - sigma_cg_inf - sigma_cq_inf )
endfunction

function g = grafico(x, y, x_label, y_label, cor)
    xgrid
    plot(x, y, cor)
    xlabel(x_label, "fontname", "times bold", "fontsize", 3)
    ylabel(y_label, "fontname", "times bold", "fontsize", 3)
endfunction

function [tabela, fck2, prot_max, prot_nec, aprov, vao] = estudo1(b, h, bf, hf, g, q, l, fck, fck_lim, aproveitamento_minimo)
    tabela   = []
    fck2     = []
    prot_max = []
    prot_nec = []
    aprov    = []
    vao      = []

    i = 1

    while fck <= fck_lim
        r = resultados(b, h, bf, hf, g, q, l, fck)

        aproveitamento = r(5)

        while aproveitamento_minimo > aproveitamento
            l = l + 0.5
            r = resultados(b, h, bf, hf, g, q, l, fck)
            aproveitamento = r(5)
        end

        tabela   = [tabela; r]
        fck2     = [fck2; r(2)]
        prot_max = [prot_max; r(3)]
        prot_nec = [prot_nec; r(4)]
        aprov    = [aprov; r(5)]
        vao      = [vao; r(6)]

        fck = fck + 5
    end
endfunction

function x = estudo2(b, h, bf, hf, g, q, l_desejado, fck)
    r = resultados(b, h, bf, hf, g, q, l_desejado, fck)

    validacao_1 = r(7)
    validacao_2 = r(8)

    while validacao_1 == 0 | validacao_2 == 0
        h = h * 1.05
        b = b * 1.05
        
        if bf > 0 then
            bf = bf * 1.05
            hf = hf * 1.05
        end

        r = resultados(b, h, bf, hf, g, q, l_desejado, fck)

        validacao_1 = r(7)
        validacao_2 = r(8)
    end
    
    x = [b h bf hf]
endfunction

//fpyk = 1710
//Ap = 9.57 * (10^(-3))
//Npd = Ap * ( fpyk / 1.15 ) * 1000 // = 14210.6

// === Primeira parte do estudo de caso ===

// Características da viga T
b = 0.2
h = 0.6
bf = 0.6
hf = 0.1
g = 15
q = 12

// 40 =< fck =< 90
l = 8
fck = 40
fck_lim = 90
aproveitamento_minimo = 90

disp("40 =< fck =< 90")
[tabela, fck2, prot_max, prot_nec, aprov, vao] = estudo1(b, h, bf, hf, g, q, l, fck, fck_lim)
disp(tabela)

// 150 =< fck =< 800
fck = 150
fck_lim = 250
aproveitamento_minimo = 94
//grafico(fck2, vao, "x", "y", 'r')

disp("150 =< fck =< 250")
[tabela, fck2, prot_max, prot_nec, aprov, vao] = estudo1(b, h, bf, hf, g, q, l, fck, fck_lim)
disp(tabela)

grafico(fck2, vao, "fck (MPa)", "l (m)", 'b')
legend(['CC e CAD';'CUAD'], ["in_lower_right"]);

// === Segunda parte do estudo de caso ===
l_desejado = 15.5
fck = 40
fck_lim = 90
est2 = []

while fck_lim >= fck

    r = estudo2(b, h, bf, hf, g, q, l_desejado, fck)
    est2 = [est2; r]

    fck = fck + 5
end

disp(est2)