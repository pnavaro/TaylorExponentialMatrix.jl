using LinearAlgebra
import LinearAlgebra:BlasFloat, checksquare

#"""
#[F,E] = log2(X) returns arrays F and E such that X=F⋅2^E
#. The values in F are typically in the range 0.5 <= abs(F) < 1.
#"""
#log2(X) = 

export expm2

"""
    expm2( A )

Matrix exponential.  ALGORITHM 1
expm2(X) is the matrix exponential of X.  expm2 is computed using
a scaling and squaring algorithm with a Taylor approximation evaluated 
using a new decomposition.
 
The choice of order and scaling is based solely on norm(A,1)
"""
function expm2(A::StridedMatrix{T}) where {T<:BlasFloat}

    n = checksquare(A)

    """
        EXPMCHK Check the class of input A and
        initialize M_VALS and THETA accordingly.
    """
    function expmchk()

        classA = eltype(A)

        if classA == Float64
            m_vals = [1 2 4 8 12 18]
            # theta_m for m=1:18.
            theta = [
                2.220446049250313e-16  # m_vals = 1
                2.580956802971767e-08  # m_vals = 2
                3.397168839976962e-04  # m_vals = 4
                4.991228871115323e-02  # m_vals = 8
                2.996158913811580e-01  # m_vals = 12
                1.090863719290036e+00  # m_vals = 18
            ]
        elseif classA == Float32
            m_vals = [1 2 4 8 12 18]
              # theta_m for m=1:7.
            theta = [
                1.192092800768788e-7   # m_vals = 1
                5.978858893805233e-04  # m_vals = 2 
                5.116619363445086e-02  # m_vals = 4
                5.800524627688768e-01  # m_vals = 8
                       #7.795113374358031e-01
                       #9.951840790004457e-01
                       #1.223479542424143e+00
                1.461661507209034  # m_vals = 12
                3.010066362817634
            ]# m_vals = 18
        else
            @error "wrong element type"
        end

        m_vals, theta, classA

    end

    """
    Improved Paterson Stockmeyer scheme, with Pade at end point
    """
    function taylor_approximant_of_degree(order)

        if (order >= 2)
            A2 = A * A
        end
        if (order > 8)
            A3 = A * A2
        end

        if order == 1

            E = I + A

        elseif order == 2

            E = I + A + A2 / 2

        elseif order == 4

            E = I + A + A2 * (I / 2 + A / 6 + A2 / 24)

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

        elseif order == 12

            # lower values of ||_||_1
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

            # Matrix products
            A6 = A3 * A3
            q31 = a01 * I + a11 * A + a21 * A2 + a31 * A3
            q61 = (b01 * I + b11 * A + b21 * A2 + b31 * A3 + b61 * A6)
            q62 = (b02 * I + b12 * A + b22 * A2 + b32 * A3 + b62 * A6)
            q63 = (b03 * I + b13 * A + b23 * A2 + b33 * A3 + b63 * A6)
            q64 = (b04 * I + b14 * A + b24 * A2 + b34 * A3 + b64 * A6)
            q91 = q31 * q64 + q63
            q18 = q61 + (q62 + q91) * q91
            E = q18

        end

        return E

    end

    m_vals, theta, classA = expmchk() # Initialization

    normA = opnorm(A, 1)

    if normA <= theta[end]
        # no scaling and squaring is required.
        s = 0
        for method_selector = 1:length(m_vals)
            if normA <= theta[method_selector]
                F = taylor_approximant_of_degree(m_vals[method_selector])
                break
            end
        end
    else
        x = normA / last(theta)
        t, s = significand(x)/2, exponent(x)+1
        s = s - Int(t == 0.5) # adjust s if normA/theta(end) is a power of 2.
        A = A ./ 2^s    # Scaling
        method_selector = length(m_vals)
        F = taylor_approximant_of_degree(m_vals[method_selector])
        for i = 1:s
            F = F * F  # Squaring
        end
    end

    @show number_of_scalings = s
    @show cost = s + method_selector - 1
    @show order = m_vals[method_selector]

    F

end
