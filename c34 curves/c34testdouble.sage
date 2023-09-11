class TestDouble(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(31^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    self.C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
  


  def test_double_0(self) :
    C_31 = self.C_31
    D1 = C_31.zero_divisor();
    D2 = C_31.zero_divisor();
    self.assertEqual(double(D1), D2)



  def test_double_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(2*D1) = 21
    D1 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 0, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, 1], [z4^3 + z4^2 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4, z4^2 + 1, 1], [z4^3 + z4 + 1, 0, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 2, 1], [1, 1, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3_3, [[2*z3, 1], [2*z3 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[2*z3, 2*z3^2 + 1, 1], [z3^2, z3, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31, [[17, 1], [24, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[3, 17, 1], [10, 3, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31_2, [[8*z2 + 22, 1], [19*z2 + 8, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[13*z2 + 8, 17*z2 + 27, 1], [15*z2 + 13, 16*z2 + 13, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_1009, [[466, 1], [801, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[356, 824, 1], [221, 932, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # Test case where type(2*D1) = 22
    D1 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[70, 1], [767, 0, 201, 0, 0, 1], []])
    self.assertEqual(2*D1, D2)




  def test_double_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(2*D1) = 41
    D1 = C34CurveDivisor(C_2, [[0, 1, 1], [1, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4^3, 1], [z4^2 + z4, z4^3 + z4, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2, z4^2, 1, 1], [z4 + 1, z4^3 + z4, z4, 0, 1], [z4^2, z4^2, z4, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3, [[0, 0, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 2, 1, 1], [2, 0, 0, 0, 1], [2, 1, 1, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3_3, [[2, z3^2 + 2*z3 + 1, 1], [2*z3 + 1, z3^2 + z3, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 1, z3, z3^2, 1], [2, 2, z3^2 + z3 + 1, 0, 1], [2*z3^2 + 2*z3 + 2, z3^2 + 1, z3^2 + 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31, [[23, 3, 1], [1, 22, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[13, 19, 18, 1], [0, 27, 24, 0, 1], [4, 5, 12, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31_2, [[6*z2 + 27, 21*z2, 1], [24*z2 + 18, 9*z2 + 18, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[5*z2 + 15, 15*z2 + 15, 22*z2 + 12, 1], [11, 27*z2 + 24, 25*z2 + 7, 0, 1], [13*z2 + 1, 11*z2 + 3, 28*z2 + 23, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_1009, [[805, 102, 1], [207, 1005, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[528, 345, 309, 1], [773, 772, 875, 0, 1], [467, 626, 890, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    # Test case where type(2*D1) = 43
    D1 = C34CurveDivisor(C_1009, [[387, 90, 1], [286, 634, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[537, 134, 1], [966, 833, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # Test case where type(2*D1) = 44
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1], [], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3_3, [[2*z3 + 2, 0, 1], [z3^2 + 2*z3 + 2, z3^2 + 2*z3, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[1], [], []])
    self.assertEqual(2*D1, D2)

    # TODO : Come up with a type 44 example over a larger curve
    #        Test case where D1 = 2*P
    C  = C34Curve(GF(997), [283, 629, 623, 533, 162, 895, 746, 521, 852])
    D1 = C34CurveDivisor(C, [[8, 553, 1], [182, 895, 0, 1], []])
    D2 = C.zero_divisor()
    self.assertEqual(2*D1, D2)



  def test_double_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(2*D1) = 43
    D1 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, z4^2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, z4^3 + z4^2 + z4, 1], [1, z4^2 + 1, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3_3, [[1, 1], [z3^2, 0, z3^2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[2*z3 + 2, 2*z3 + 2, 1], [2*z3, 2, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31, [[16, 1], [12, 0, 30, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[18, 1, 1], [22, 19, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31_2, [[25*z2 + 1, 1], [17*z2 + 12, 0, 15*z2 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[7*z2 + 6, 4*z2 + 16, 1], [18*z2 + 26, 13*z2 + 22, 0, 1], []])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_1009, [[714, 1], [692, 0, 574, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[271, 641, 1], [528, 67, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # Test case where type(2*D1) = 42
    C = C34Curve(GF(997), [278, 360, 340, 309, 34, 816, 532, 458, 688])
    D1 = C34CurveDivisor(C, [[487, 1], [982, 0, 748, 0, 0, 1], []])
    D2 = C34CurveDivisor(C, [[487, 1], [744, 0, 1], []])



  def test_double_31(self) :
    C_2, C_3, C_11, C_31, C_1009 = self.C_2, self.C_3, self.C_11, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where D is typical and type(2D) = 61
    D1 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^2, z4^3 + z4 + 1, 1], [z4 + 1, z4^2 + 1, z4^3, 0, 1], [z4^2 + 1, z4^3 + z4^2, z4^2 + z4 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[0, z4^3 + z4^2 + 1, z4, 1], [z4^2 + z4, z4^3 + z4^2, z4^3 + z4^2 + 1, 0, 1], [0, z4^3 + z4 + 1, z4^3 + z4^2, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [2, 1, 2, 0, 1], [0, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[0, 2, 2, 1], [1, 2, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 2*z3^2 + 1, 2*z3, 1], [z3^2 + z3, 2*z3 + 2, 2*z3^2 + z3 + 2, 0, 1], [z3^2 + 2*z3 + 1, z3, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2, 2, z3 + 2, 1], [0, z3^2 + 2*z3, z3^2 + z3, 0, 1], [z3 + 2, z3^2 + z3 + 2, 2*z3, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31, [[3, 14, 5, 1], [27, 23, 30, 0, 1], [1, 21, 8, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[26, 20, 3, 1], [7, 7, 25, 0, 1], [0, 25, 16, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_31_2, [[5*z2 + 15, 26*z2 + 12, 10*z2, 1], [25*z2 + 5, 22*z2 + 10, 28*z2 + 27, 0, 1], [24*z2 + 14, 29*z2 + 23, 26*z2 + 27, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[23*z2 + 13, 25*z2 + 18, 24*z2 + 27, 1], [17*z2 + 12, 3, 10*z2 + 29, 0, 1], [21*z2 + 8, 4*z2 + 10, 4*z2 + 30, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    D1 = C34CurveDivisor(C_1009, [[109, 761, 961, 1], [384, 439, 627, 0, 1], [979, 769, 64, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[491, 684, 797, 1], [432, 887, 687, 0, 1], [305, 959, 894, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    # D is typical and type(2D) = 62
    # TODO : Find an example over a larger field
    C = C34Curve(GF(97), [13, 17, 92, 0, 5, 80, 69, 22, 60])
    D1 = C34CurveDivisor(C, [[38, 33, 83, 1], [59, 83, 47, 0, 1], [76, 72, 61, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[10, 1], [52, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # D is typical and type(2D) = 63
    D1 = C34CurveDivisor(C_1009, [[896, 14, 20, 1], [577, 476, 537, 0, 1], [947, 896, 286, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[335, 875, 1], [616, 976, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # D is typical and type(2D) = 64
    # TODO : Find an example over a larger field
    C = C34Curve(GF(7), [4, 4, 4, 5, 3, 0, 0, 6, 1])
    D1 = C34CurveDivisor(C, [[1, 0, 4, 1], [4, 0, 4, 0, 1], [4, 6, 6, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[5, 1], [2, 0, 6, 0, 0, 1], []])
    self.assertEqual(2*D1, D2)

    # D is typical and type(2D) = 65
    # TODO : Find an example over a larger field
    C = C34Curve(GF(37), [30, 14, 36, 25, 12, 30, 6, 32, 2])
    D1 = C34CurveDivisor(C, [[12, 0, 28, 1], [16, 23, 9, 0, 1], [15, 24, 25, 0, 0, 1]])
    D2 = C.zero_divisor()
    self.assertEqual(2*D1, D2)

    # D is atypical and type(2D - G) = 51
    C = C34Curve(GF(997), [221, 871, 509, 149, 696, 991, 69, 850, 347])
    D1 = C34CurveDivisor(C, [[293, 788, 0, 1], [147, 851, 238, 0, 1], [779, 61, 17, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[84, 307, 49, 1], [915, 967, 830, 0, 1], [50, 226, 814, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    # D is atypical and type(2D - G) = 52
    # This case is impossible, since D = G + flip(P) for some point P
    # If 2D - G is type 52, then 2D - G is equivalent to some point Q and
    #
    #   2D - G == Q   ==>   D - P == Q   ==>   D = Q + P

    # D is atypical and type(2D - G) = 53
    C = C34Curve(GF(997), [256, 805, 462, 853, 250, 490, 487, 360, 639])
    D1 = C34CurveDivisor(C, [[814, 450, 0, 1], [353, 930, 709, 0, 1], [79, 158, 199, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[875, 153, 299, 1], [893, 66, 38, 0, 1], [800, 799, 121, 0, 0, 1]])
    self.assertEqual(2*D1, D2)

    # D is atypical and type(2D - G) = 54
    # This case is also impossible.
    # D = G + flip(P) for some point P.
    # If 2D - G is type 54, then 2D - G is equivalent to flip(Q) for some point Q
    #   flip(Q) = 2D - G = D + flip(P)
    #   D = P + flip(Q)
    # Since the decomposition of D into the sum of type 11 and 22 divisors is unique
    #   G = P and flip(Q) = flip(P)   ==>   D = G + flip(P) = P + flip(P) = 0

