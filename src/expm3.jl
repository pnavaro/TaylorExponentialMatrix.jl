using LinearAlgebra
import ExpmV.normAm


export expm3

"""
    expm3(X)   

Matrix exponential computed using
a scaling and squaring algorithm with a Taylor approximation evaluated 
using a new decomposition.
   
When the maximal order is selected, estimates of norm(A^9,1) are used
to reduce overscaling

The choice of order and scaling is based solely on norm(A,1)
"""
function expm3(A)

    s, degree, Tpowers = expm_params(A)

    @show s
    # Rescale the powers of A appropriately.

    if s != 0
        for (k, v) in Tpowers
            @show k
            Tpowers[k] ./= 2 .^ (s * k)
        end
    end

    # Evaluate the Taylor polynomial to given degree.
    E  = zero(A)
    E .= taylor_approximant_of_degree(degree, Tpowers)

    for i = 1:s
        E .= E * E  # Squaring
    end

    E

end

"""
    expm_params(A)

expm_helper Obtain scaling parameter and order of the Taylor 
approximant
"""
function expm_params(A)

    s :: Int = 0

    theta = getThetaSet(A)
    @show normA = opnorm(A, 1)

    Tp = Dict()

    # Find appropriate order:
    Tp[1] = A
    ## Check if degree 1 is sufficient (No info available)
    if (normA < theta[1])
        degree = 1
        return s, degree, Tp
    end
    Tp[2] = A * A
    if (normA < theta[2])
        degree = 2
        return s, degree, Tp
    end
    if (normA < theta[3])
        degree = 4
        return s, degree, Tp
    end
    # At degree 8, we could do the first estimates
    # with pairs (2,3), (3,4), (3,5), (2,9) [Also (2,5), (2,7)]

    if (normA < theta[4])
        degree = 8
        return s, degree, Tp
    end
    Tp[3] = A * Tp[2]
    if (normA < theta[5])
        degree = 12
        return s, degree, Tp
    end

    # Otherwise degree 18
    # Here, we only gain the pairs (3, 8) and (4, 7), 
    # but both are not expected to
    # yield improvements, however, we also get (2, 2k+1) until (2,17)
    degree = 18
    Tp[6] = Tp[3] * Tp[3]
    d2 = norm(Tp[2], 1)^(1 / 2)
    d3 = norm(Tp[3], 1)^(1 / 3)
    d6 = norm(Tp[6], 1)^(1 / 6)
    eta = min(max(d2, d3), normA)
    quotients = minimum([d2 / normA, d3 / normA, d6 / normA])
    if (quotients < 1 / 2^4)  ## strong decay
        d9 = normAm(A, 9)^(1 / 9)
        d19 = normA
        eta = min([eta, max(d2, d9), max(d2, d19)])
    end
    s = max(ceil(log2(eta / theta[6])), 0)
    if isinf(s)
        # Overflow in ell subfunction. Revert to old estimate.
        x = opnorm(T, 1) / theta[end]
        t, s = significand(x) / 2, exponent(x) + 1
        s = s - Int(t == 0.5) # adjust s if normA/theta(end) is a power of 2.
    end

    return s, degree, Tp
end

"""
    taylor_approximant_of_degree(order, Tp)

Improved Paterson Stockmeyer scheme, with Pade at end point
"""
function taylor_approximant_of_degree(order, Tp)

    A = Tp[1]
    if order >= 4
        A2 = Tp[2]
    end
    if order > 8
        A3 = Tp[3]
    end

    if order == 1
        E = I + A
    elseif order == 2
        E = I + A
    elseif order == 4
        E = I + A + A2 / 2 * (I + A / 3 + A2 / 3 / 4)
    elseif order == 8
        # Minimizes ||coefficients||_1
        x3 = 2 / 3
        a1 = 1 / 88 * (1 + sqrt(177)) * x3
        a2 = 1 / 352 * (1 + sqrt(177)) * x3
        u0 = 1
        u1 = 1
        u2 = 1 / 630 * (857 - 58 * sqrt(177))
        c0 = (-271 + 29 * sqrt(177)) / (315 * x3)
        c1 = (11 * (-1 + sqrt(177))) / (1260 * x3)
        c2 = (11 * (-9 + sqrt(177))) / (5040 * x3)
        c4 = -((-89 + sqrt(177)) / (5040 * x3^2))
        # Matrix products
        A4 = A2 * (a1 * A + a2 * A2)
        A8 = (x3 * A2 + 1 * A4) * (c0 * I + c1 * A + c2 * A2 + c4 * A4)
        E = u0 * I + u1 * A + u2 * A2 + +A8
    elseif order == 12 # lower values of ||_||_1
        a01 = -0.0186023205146205532243437300433
        a02 = 4.60000000000000000000000000000
        a03 = 0.211693118299809442949323323336
        a04 = 0
        a11 = -0.00500702322573317730979741843919
        a12 = 0.992875103538486836140479571505
        a13 = 0.158224384715726725371768893252
        a14 = -0.131810610138301840156819349464
        a21 = -0.573420122960522263905952420789
        a22 = -0.132445561052799638845074997454
        a23 = 0.165635169436727415011171668419
        a24 = -0.0202785554058925907933568229945
        a31 = -0.133399693943892059700768926983
        a32 = 0.00172990000000000000000000000000
        a33 = 0.0107862779315792425026320640108
        a34 = -0.00675951846863086359778560766482
        q31 = a01 * I + a11 * A + a21 * A2 + a31 * A3
        q32 = a02 * I + a12 * A + a22 * A2 + a32 * A3
        q33 = a03 * I + a13 * A + a23 * A2 + a33 * A3
        q34 = a04 * I + a14 * A + a24 * A2 + a34 * A3
        # Matrix products
        q61 = q33 + q34 * q34
        E = (q31 + (q32 + q61) * q61)
    elseif order == 18
        # Minimizes ||coefficients||_1
        a01 = 0
        a11 = -0.100365581030144620014614939405
        a21 = -0.00802924648241156960116919515244
        a31 = -0.000892138498045729955685466128049
        b01 = 0
        b11 = 0.3978497494996450761451961277102845756965081084076856223845951607640145373149032030404660339703426170
        b21 = 1.367837784604117199225237068782228242106453540654915795267462332707000476284638745738812082761860458
        b31 = 0.4982896225253826775568588172622733561849181397319696269923450773179744578906675427707618377504305561
        b61 = -0.0006378981945947233092415500564919285518773827581013332055148653241353120789646323186965398317523194760
        b02 = -10.96763960529620625935174626753684124863041876254774214673058385106461743913502064396554589372626845
        b12 = 1.680158138789061971827854017248105095278579593010566858091585875627364747895724070033586802947436157
        b22 = 0.05717798464788655127028717132252481317274975093707600394506105236273081373356534970690710643085727120
        b32 = -0.006982101224880520842904664838015531664389838152077753258353447125605205335288966277257918925881337834
        b62 = 0.00003349750170860705383133673406684398020225996444991565389728295589367037178816169580298011831485225359
        b03 = -0.09043168323908105619714688770183495499470866281162649735086602288456671216199491949073419844120202066
        b13 = -0.06764045190713819075600799003797074854562826013273060070581796689746210078264908587143648477465431597
        b23 = 0.06759613017704596460827991950789370565071000768138331455747136830953224978586305376886338214283464385
        b33 = 0.02955525704293155274260691822422312437473046472743510159951445245685893720550532525991666174704105350
        b63 = -0.00001391802575160607011247399793437424028074102305234870853220069004278819595344896368546425681168813708
        b04 = 0
        b14 = 0
        b24 = -0.0923364619367118592764570775143
        b34 = -0.0169364939002081717191385115723
        b64 = -0.0000140086798182036159794363205726

        # Extra matrix product
        A6 = Tp[6]
        q31 = a01 * I + a11 * A + a21 * A2 + a31 * A3
        q61 = (b01 * I + b11 * A + b21 * A2 + b31 * A3 + b61 * A6)
        q62 = (b02 * I + b12 * A + b22 * A2 + b32 * A3 + b62 * A6)
        q63 = (b03 * I + b13 * A + b23 * A2 + b33 * A3 + b63 * A6)
        q64 = (b04 * I + b14 * A + b24 * A2 + b34 * A3 + b64 * A6)
        q91 = q31 * q64 + q63
        q18 = q61 + (q62 + q91) * q91
        E = q18
    else
        throw("expm3= degree:$order")
    end

    return E

end

function getThetaSet(A)

    classA = eltype(A)
    m_vals = [1 2 4 8 12 18]
    if classA == Float64

        # theta_m for m=2:18.
        theta = [
            2.220446049250313e-16  # m_vals = 1
            2.580956802971767e-08  # m_vals = 2
            3.397168839976962e-04  # m_vals = 4
            4.991228871115323e-02  # m_vals = 8
            2.996158913811580e-01  # m_vals = 12
            1.090863719290036e+00
        ]# m_vals = 18

    elseif classA == Float32

        # theta_m for m=1:7.
        theta = [
            1.192092800768788e-7   # m_vals = 1
            5.978858893805233e-04  # m_vals = 2 
            5.116619363445086e-02  # m_vals = 4
            5.800524627688768e-01  # m_vals = 8
           #7.795113374358031e-01
           #9.951840790004457e-01
           #1.223479542424143e+00
            1.461661507209034e+00  # m_vals = 12
            3.010066362817634e+00
        ]# m_vals = 18

    else
        throw("Input matrix element type must be Float32 or Float64")
    end

    return theta

end
