function x = resultados(b, h, bf, hf, g, q, l, fck, Npd)
    compressao_maxima = 0.6 * fck * 1000

    [Ac, I, Wsup, Winf] = caractGeometricas(b, h, bf, hf)
    gg = g + (Ac * 25)
    [Mgmax, Mqmax, sigma_cg_inf, sigma_cg_sup, sigma_cq_inf, sigma_cq_sup] = esforcosSolicitantes(gg, q, l, Wsup, Winf)

    //Npd = Ap * ( fpyk / 1.15 ) * 1000

    sigma_c = Npd / 0.16

    Np = ( - sigma_cg_inf - sigma_cq_inf ) * Ac

    sigma_c_min = - Npd + sigma_cg_sup + sigma_cq_sup

    if abs(sigma_c_min) < abs(compressao_maxima) then
        v1 = 1
    else
        v1 = 0
    end

    if abs(Npd) < abs(compressao_maxima) then
        v2 = 1
    else
        v2 = 0
    end

//    v3 = v1 * v2

    if (v1 * v2) == 1 then
        x = [Npd fck sigma_c_min abs(compressao_maxima) (abs(sigma_c_min)) ((abs(sigma_c_min) / abs(compressao_maxima)) * 100) l fpyk Ap]
    end
endfunction

function [Ac, I, Wsup, Winf] = caractGeometricas(b, h, bf, hf)
// Seção retangular
    if bf == 0 & hf == 0 then
        Ac = b * h
        I  = (b * (h^3)) / 12
        ysup = h / 2
        yinf = h / 2
        Wsup = I / ysup
        Winf = I / yinf
    else
        h = (h - hf)
        A1 = bf * hf
        A2 = b  * h
        Ac = A1 + A2
        Msx = (A1 * (h + (hf / 2))) + (A2 * (h / 2))
        yg = Msx / Ac
        I_xg1 = ((bf * (hf^3)) / 12) + (bf * hf * (((h + (hf / 2)) - yg)^2))
        I_xg2 = ((b * (h^3)) / 12) + (b * h * ((yg - (h / 2))^2))
        I_xg = I_xg1 + I_xg2
        I  = I_xg
        ysup = (h + hf) - yg
        yinf = yg
        Wsup = I / ysup
        Winf = I / yinf
    end
    
//    disp("Área da seção (Ac): " + Ac + " m^2")
//    disp("Inércia da seção (Ic): " + I + " m^4")
//    disp("Módulo resistente superior (Wc_sup): " + Wsup + " m^3")
//    disp("Módulo resistente inferior (Wc_inf): " + Winf + " m^3")  
endfunction

function [Mgmax, Mqmax, sigma_cg_inf, sigma_cg_sup, sigma_cq_inf, sigma_cq_sup] = esforcosSolicitantes(g, q, l, Wsup, Winf)
    j = (l^2 / 8)
    Mgmax = g * j
    Mqmax = q * j
    
    sigma_cg_inf = + Mgmax / Winf
    sigma_cg_sup = - Mgmax / Wsup
    sigma_cq_inf = + Mqmax / Winf
    sigma_cq_sup = - Mqmax / Wsup
endfunction

b = 0.2
h = 0.6
bf = 0.6
hf = 0.1
l = 16
g = 11
q = 12
fck = 40
Npd = 14210.6
fpyk = 1710
Ap = 9.57 * (10^(-3))

x = []
eixox = []
vao = []
protensao = []

i = fck
ii = 1

while i <= 250
    x = resultados(b, h, bf, hf, g, q, l, i, Npd)
    if x(6) < 92 then
        l = l + 1
        x = resultados(b, h, bf, hf, g, q, l, i, Npd)
    end

    if x(2) == (i-5) then
        x = resultados(b, h, bf, hf, g, q, l, i, Npd)
    end

    eixox(ii) = x(2)
    vao(ii) = x(7) 
    protensao(ii) = x(5) / 1000

    disp([x(2) x(4) x(5) x(6) x(7)])
    i = i + 5
    ii = ii + 1
end

//plot(eixox, vao)
xgrid
//xlabel("$\Large{fck} \quad \text{(MPa)}$", "fontname", "times bold", "fontsize", 3)
//ylabel("$\Large{l} \quad \text{(m)}$", "fontname", "times bold", "fontsize", 3)

plot(vao, protensao)
xlabel("$\Large{l} \quad \text{(m)}$", "fontname", "times bold", "fontsize", 3)
ylabel("$\Large{\sigma_{c,min}} \quad \text{(MPa)}$", "fontname", "times bold", "fontsize", 3)
