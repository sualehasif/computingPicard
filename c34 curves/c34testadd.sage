load("c34crv.sage")

class TestAdd(unittest.TestCase) :
  def setUp(self) :
    self.z2 = GF(31^2).gen()
    self.z3 = GF(3^3).gen()
    self.z4 = GF(2^4).gen()
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    self.C_2 = C34Curve(GF(2), [0, 0, 1, 0, 1, 1, 1, 1, 1])
    self.C_3 = C34Curve(GF(3), [0, 2, 1, 1, 0, 2, 2, 0, 2])
    self.C_11 = C34Curve(GF(11), [6, 5, 1, 8, 6, 3, 7, 10, 4]) # The vertical tangent line at P = (9 : 6 : 1) intersects P with multiplicity 2
    self.C_31 = C34Curve(GF(31), [2, 29, 13, 7, 28, 25, 17, 13, 17])
    #self.C_41 = C34Curve(GF(41), [26, 22, 16, 37, 30, 22, 7, 22, 29]) # Vertical tangents at (31 : 33 : 1) and (28 : 22 : 1) intersect with m = 2
    self.C_2_4 = C34Curve(GF(2^4), [z4^3 + 1, z4, z4^3 + z4^2 + 1, z4^2, z4^3 + z4, z4^3 + z4 + 1, z4^3 + z4, z4^3 + z4^2 + z4, z4^3 + 1])
    self.C_3_3 = C34Curve(GF(3^3), [z3 + 2, z3 + 1, z3^2 + z3 + 2, 2*z3^2 + z3 + 1, 2*z3^2 + z3, z3^2 + 2, 2*z3^2 + z3, 2*z3^2 + 1, 2])
    self.C_31_2 = C34Curve(GF(31^2), [9*z2 + 11, 11*z2 + 28, 17*z2 + 2, 24*z2 + 22, 6*z2 + 29, 21*z2 + 8, 18*z2 + 5, 18*z2 + 15, 13*z2 + 10])
    self.C_1009 = C34Curve(GF(1009), [715, 286, 748, 7, 92, 446, 122, 935, 314])
  


  def test_add_11_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case, type(lcm(D1, D2)) = 21, over several fields
    D1 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^2 + z4, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2, 1],
                                 [z4^2 + z4 + 1, z4 + 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 1, 1], [0, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 1, 1], [z3 + 2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3 + 2, 0, 1], [z3^2 + 2*z3, 2, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[22, 1], [13, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[28, 1], [17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 20, 1], [27, 19, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[17*z2 + 25, 1], [13*z2 + 9, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[16*z2 + 4, 1], [14*z2 + 18, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[13*z2 + 4, 5*z2 + 10, 1], [20*z2 + 28, 2*z2 + 29, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[202, 1], [648, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[861, 1], [413, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[743, 475, 1], [374, 54, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    self.assertEqual(add_11_11(D1, D2), D3)

    # Test rarer case, type(lcm(D1, D2)) = 22
    D1 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + 1, 1], [z4^2 + 1, 0, z4^2, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [z3 + 1, 0, 2*z3 + 1, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[11, 1], [26, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[11, 1], [28, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[11, 1], [15, 0, 23, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [26*z2 + 7, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [30*z2 + 6, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[27*z2 + 9, 1], [4*z2 + 27, 0, 25*z2 + 13, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[70, 1], [767, 0, 201, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)



  def test_add_21_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case, type(lcm(D1, D2)) = 31, over several fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4^2, 1],
                                 [z4^2 + z4 + 1, z4 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^2 + z4, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4 + 1, 1, z4^2 + z4 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^3 + z4, z4^3 + 1, 0, 1],
                                 [z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[2, 0, 1], [2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[2, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[2, 0, 0, 1], [1, 2, 2, 0, 1], [1, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3 + 1, 2*z3^2, 1], [2*z3^2 + 1, 2*z3^2 + z3, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, 1], [2*z3^2 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, 2*z3^2 + 2, z3^2 + 2*z3 + 2, 1],
                                 [z3^2 + z3, 2*z3^2 + 2*z3, z3^2 + 2, 0, 1],
                                 [2*z3, z3^2 + 2*z3 + 1, 2*z3^2 + z3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[1, 10, 1], [14, 24, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[30, 1], [25, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[19, 12, 5, 1], [26, 26, 11, 0, 1], [13, 21, 9, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[22*z2 + 11, 18*z2 + 9, 1], [29*z2 + 2, 14*z2 + 29, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[28*z2 + 30, 1], [5*z2 + 24, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[16*z2 + 28, 9*z2 + 8, 28*z2 + 21, 1], [15*z2 + 24, 2*z2 + 20, 2*z2 + 20, 0, 1], [6*z2 + 1, 21*z2 + 29, 11*z2 + 25, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_1009, [[118, 576, 1], [929, 39, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[958, 1], [670, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[930, 232, 419, 1],
                                  [135, 568, 765, 0, 1],
                                  [291, 230, 72, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 32
    D1 = C34CurveDivisor(C_1009, [[255, 627, 1], [818, 275, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[512, 1], [93, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[519, 1], [203, 0, 189, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where D2 < D1, but D1 =/= 2*D2
    D1 = C34CurveDivisor(C_1009, [[247, 321, 1], [312, 320, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[25, 1], [294, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[341, 709, 523, 1], [641, 645, 645, 0, 1], [45, 337, 341, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = 2*D2
    D1 = C34CurveDivisor(C_1009, [[124, 767, 1], [296, 225, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[617, 1], [106, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[741, 512, 28, 1], [553, 948, 330, 0, 1], [665, 955, 379, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

  
  
  def test_add_21_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[0, z4^2 + z4, 1], [z4^3 + z4^2 + z4, z4^2, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + 1, z4^3 + z4^2 + z4 + 1, 1], [z4 + 1, z4^2, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4 + 1, z4^3 + z4^2 + 1, 1],
                                 [0, z4^2 + 1, z4^3 + 1, 0, 1],
                                 [z4^3 + z4, z4^3 + z4, z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[1, 2, 1], [0, 1, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 0, 1, 1], [0, 0, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3, z3^2 + 2*z3 + 1, 1], [z3 + 1, z3^2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 2*z3^2 + z3 + 1, 1],
                                 [1, z3^2 + z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3, 2*z3^2 + 2*z3 + 1, 2*z3 + 2, 1],
                                 [2*z3, 2*z3^2 + 1, z3^2, 0, 1],
                                 [z3, z3^2 + 2*z3 + 1, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[10, 17, 1], [15, 14, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[29, 25, 1], [8, 17, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[30, 13, 29, 1], [18, 11, 13, 0, 1], [21, 15, 27, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[9*z2 + 19, 17*z2 + 9, 1], [14*z2 + 22, 28*z2 + 18, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[27, 13*z2 + 27, 1], [6*z2 + 24, 25*z2 + 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[6*z2 + 13, 3*z2 + 13, 11*z2 + 6, 1],
                                  [20*z2 + 15, 28*z2 + 17, 3*z2 + 1, 0, 1],
                                  [16*z2 + 16, 5*z2 + 22, 19*z2 + 9, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[688, 179, 1], [360, 620, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[13, 218, 1], [131, 369, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[754, 638, 275, 1],
                                  [904, 369, 1007, 0, 1],
                                  [104, 959, 281, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    
    # Test case where type(lcm(D1, D2)) = 43
    D1 = C34CurveDivisor(C_1009, [[867, 159, 1], [511, 196, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[608, 796, 1], [840, 614, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[499, 368, 1], [229, 871, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    # Test case where type(lcm(D1, D2)) = 44
    D1 = C34CurveDivisor(C_1009, [[103, 411, 1], [993, 797, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[103, 411, 1], [203, 292, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 31
    D1 = C34CurveDivisor(C_1009, [[903, 416, 1], [788, 201, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[882, 42, 1], [617, 615, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[462, 653, 489, 1], [61, 262, 212, 0, 1], [143, 489, 926, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 32
    D1 = C34CurveDivisor(C_31, [[18, 21, 1], [26, 22, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[18, 21, 1], [4, 29, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[10, 13, 0, 1], [23, 22, 25, 0, 1], [15, 28, 14, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)



  def test_add_21_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1, 1], [1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^2 + z4, 1],
                                 [z4^3 + z4 + 1, z4^2 + z4 + 1, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, 1], [z4^3, 0, z4^3 + z4^2 + z4, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4, z4 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^3 + z4^2 + z4, z4^2, 0, 1],
                                 [z4^3 + z4^2 + 1, 1, z4^2 + z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 0, 1], [2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 2, 1], [2, 1, 1, 0, 1], [2, 1, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3, 2*z3^2 + 2*z3 + 1, 1],
                                 [z3^2, z3^2 + 2, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 1, 1], [z3^2 + 2*z3, 0, 2*z3 + 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3 + 1, 2*z3 + 1, 1],
                                 [z3, z3 + 1, 2*z3^2 + z3 + 1, 0, 1],
                                 [z3^2 + z3, 2*z3^2 + 2*z3, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[7, 23, 1], [8, 5, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[15, 1], [6, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[0, 18, 3, 1], [15, 20, 21, 0, 1], [15, 11, 10, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[28*z2 + 9, 24*z2 + 5, 1], [27*z2 + 3, 29*z2 + 3, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[16*z2 + 14, 1], [7*z2 + 15, 0, 29*z2 + 27, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[z2 + 29, 29*z2 + 3, 8*z2 + 13, 1],
                                  [3*z2 + 10, 2*z2 + 20, 13*z2 + 1, 0, 1],
                                  [18*z2 + 13, 14*z2 + 23, 3*z2 + 29, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[35, 881, 1], [1006, 1007, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[858, 1], [369, 0, 444, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[809, 604, 751, 1],
                                  [535, 20, 774, 0, 1],
                                  [101, 384, 817, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    # Test case where type(lcm(D1, D2)) = 42
    D1 = C34CurveDivisor(C_1009, [[811, 316, 1], [368, 408, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[90, 1], [715, 0, 824, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[318, 1], [214, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 31
    D1 = C34CurveDivisor(C_1009, [[441, 228, 1], [722, 896, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[94, 1], [333, 0, 429, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[662, 853, 899, 1], [118, 517, 129, 0, 1], [307, 1000, 296, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)



  def test_add_22_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test case where type(lcm(D1, D2)) = 31 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 0, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 0, 0, 1], [1, 1, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^2 + z4, 1], [z4^3 + z4, 0, z4^2 + z4, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2, z4^3 + 1, 0, 1],
                                 [z4^3 + z4^2 + z4 + 1, z4^3 + z4 + 1, z4^2 + z4, 0, 1],
                                 [z4^2 + 1, z4^3 + z4 + 1, z4^2 + z4, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 1], [0, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [0, 2, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [z3^2, 0, 2*z3^2 + 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 2, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3, 2*z3^2 + z3 + 1, 0, 1],
                                 [2*z3^2 + z3 + 1, 2, z3^2 + 2*z3 + 2, 0, 1],
                                 [2*z3^2 + z3, z3 + 2, 2*z3^2 + 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[30, 1], [1, 0, 17, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[19, 1], [18, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[12, 18, 0, 1], [13, 18, 30, 0, 1], [14, 18, 17, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[5*z2 + 18, 1], [22*z2 + 30, 0, 30*z2 + 3, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 6, 1], [9*z2 + 25, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[20*z2 + 27, 29*z2 + 24, 0, 1],
                                  [5*z2 + 5, 9*z2 + 25, 5*z2 + 18, 0, 1],
                                  [17*z2 + 6, 28*z2 + 22, 30*z2 + 3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[341, 1], [435, 0, 354, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[566, 1], [794, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[287, 907, 0, 1],
                                  [342, 794, 341, 0, 1],
                                  [485, 299, 354, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test case where type(lcm(D1, D2)) = 33
    D1 = C34CurveDivisor(C_1009, [[246, 1], [124, 0, 63, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[246, 1], [832, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[1], [], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where D2 < D1, but D1 =/= 2*D2
    D1 = C34CurveDivisor(C_1009, [[589, 1], [710, 0, 862, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[589, 1], [1004, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[834, 169, 0, 1], [82, 1004, 589, 0, 1], [162, 501, 862, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = 2*D2
    D1 = C34CurveDivisor(C_1009, [[70, 1], [767, 0, 201, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[70, 1], [605, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[864, 140, 0, 1], [981, 605, 70, 0, 1], [262, 209, 201, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)




  def test_add_22_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 43 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2, [[0, 1], [1, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 0, 1], [0, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4 + 1, 1], [z4^2, 0, z4^3 + z4 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_2_4, [[z4^2 + z4 + 1, 1], [z4^3 + z4 + 1, 0, z4^2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^3 + 1, z4^3 + z4^2 + z4, 1],
                                 [z4^3 + z4^2 + 1, 1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3, [[0, 1], [1, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 1], [0, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_3_3, [[1, 1], [0, 0, 2*z3 + 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2 + 2*z3 + 1, z3, 1],
                                 [z3^2 + 2*z3 + 1, z3^2 + 2*z3 + 2, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[19, 1], [23, 0, 21, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31, [[25, 1], [28, 0, 14, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 21, 1], [30, 13, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[17*z2 + 18, 1], [21*z2 + 14, 0, 10*z2 + 9, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_31_2, [[21*z2 + 1, 1], [3*z2 + 27, 0, 21*z2 + 26, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[27*z2 + 10, 27*z2 + 21, 1],
                                  [26*z2 + 9, 11*z2 + 22, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[67, 1], [81, 0, 579, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[697, 1], [870, 0, 233, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[836, 464, 1], [424, 51, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 33
    D1 = C34CurveDivisor(C_1009, [[350, 1], [356, 0, 171, 0, 0, 1], []])
    D2 = C34CurveDivisor(C_1009, [[350, 1], [148, 0, 107, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[350, 1], [760, 0, 1], []])
    self.assertEqual(D1 + D2, D3)



  def test_add_31_11(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(lcm(D1, D2)) = 41 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4, z4 + 1, z4^2 + 1, 1],
                                 [z4^3 + z4^2 + z4, z4^2 + 1, z4^2 + z4, 0, 1],
                                 [z4^2, z4^3 + z4^2 + z4, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, 1], [1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[1, z4^3 + z4^2 + z4 + 1, 1, 1],
                                 [z4^3 + z4^2 + 1, z4^3, z4^3 + z4 + 1, 0, 1],
                                 [z4^3 + 1, z4^3 + z4, z4 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 2, 2, 1], [2, 1, 1, 0, 1], [2, 1, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1], [2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3 + 1, 2*z3^2 + 2*z3, z3^2 + 2, 1],
                                 [2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3, z3^2 + 2, 0, 1],
                                 [z3, z3^2 + 1, 2*z3 + 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, 1], [2*z3 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2, 2*z3^2 + 1, 2*z3^2 + z3 + 2, 1],
                                 [z3^2 + z3, 2, z3^2 + z3, 0, 1],
                                 [z3^2 + 1, 2*z3^2 + z3 + 2, z3^2 + 2*z3 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[15, 4, 17, 1], [23, 4, 29, 0, 1], [10, 0, 22, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[15, 1], [16, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[18, 13, 7, 1], [29, 4, 14, 0, 1], [10, 26, 13, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[26*z2 + 15, 2*z2 + 10, 9*z2 + 26, 1],
                                  [11*z2 + 17, 27*z2 + 14, 14*z2 + 11, 0, 1],
                                  [z2 + 28, 28*z2 + 11, 16*z2 + 14, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[2*z2 + 30, 1], [14*z2 + 28, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[4*z2 + 24, 28*z2 + 13, 20*z2 + 8, 1],
                                  [10, 24*z2 + 30, 29*z2 + 7, 0, 1],
                                  [18*z2 + 11, 12*z2 + 9, 10*z2 + 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[221, 507, 55, 1],
                                  [112, 772, 678, 0, 1],
                                  [38, 507, 389, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[448, 1], [82, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[296, 653, 266, 1],
                                  [802, 145, 938, 0, 1],
                                  [482, 789, 878, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 42
    D1 = C34CurveDivisor(C_31, [[14, 8, 0, 1], [12, 1, 12, 0, 1], [6, 14, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[12, 1], [6, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[27, 1], [1, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 43
    D1 = C34CurveDivisor(C_1009, [[198, 345, 749, 1], [183, 812, 151, 0, 1], [240, 542, 986, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[75, 1], [667, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[670, 427, 1], [715, 669, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    # Test the many cases where D2 < D1
    # Test case where D1 = D2 + A, lcm(D2, A) = 0, type(A) = 21, D1 typical
    D1 = C34CurveDivisor(C_1009, [[520, 815, 68, 1], [207, 760, 607, 0, 1], [371, 0, 98, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[607, 1], [347, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[113, 837, 665, 1], [185, 516, 210, 0, 1], [593, 762, 391, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = D2 + A, lcm(D2, A) = 0, type(A) = 21, D1 atypical, D1.f has distinct roots
    D1 = C34CurveDivisor(C_1009, [[667, 301, 0, 1], [449, 18, 81, 0, 1], [476, 295, 587, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[81, 1], [740, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[978, 844, 390, 1], [895, 205, 799, 0, 1], [372, 961, 744, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = D2 + A, lcm(D2, A) = 0, type(A) = 22, D1 atypical, D1.f has distinct roots
    D1 = C34CurveDivisor(C_1009, [[478, 628, 0, 1], [375, 248, 624, 0, 1], [871, 353, 494, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[4, 1], [248, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[364, 399, 36, 1], [1007, 519, 203, 0, 1], [662, 264, 545, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = D2 + A, lcm(D2, A) = 0, type(A) = 21, D1 atypical, D1.f has a double root
    D1 = C34CurveDivisor(C_1009, [[635, 391, 0, 1], [794, 915, 700, 0, 1], [439, 358, 238, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[700, 1], [332, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[505, 411, 1], [507, 698, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = D2 + A, lcm(D2, A) = 0, type(A) = 22, D1 atypical, D1.f has a double root
    D1 = C34CurveDivisor(C_1009, [[635, 391, 0, 1], [794, 915, 700, 0, 1], [439, 358, 238, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[700, 1], [915, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[324, 578, 788, 1], [978, 376, 209, 0, 1], [406, 982, 1002, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = 2*D2 + A, lcm(D2, A) = 0, type(2*D2) = 21 
    D1 = C34CurveDivisor(C_1009, [[480, 106, 660, 1], [231, 564, 989, 0, 1], [999, 776, 752, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[596, 1], [963, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[170, 2, 86, 1], [516, 646, 230, 0, 1], [275, 830, 413, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = 2*D2 + A, lcm(D2, A) = 0, type(2*D2) = 22 
    D1 = C34CurveDivisor(C_1009, [[78, 775, 0, 1], [231, 1000, 647, 0, 1], [780, 507, 52, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[647, 1], [26, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[729, 669, 821, 1], [233, 563, 951, 0, 1], [328, 717, 308, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D1 = 3*D2
    D1 = C34CurveDivisor(C_1009, [[520, 109, 77, 1], [75, 287, 237, 0, 1], [189, 581, 596, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[788, 1], [967, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[937, 690, 568, 1], [768, 343, 384, 0, 1], [836, 60, 136, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)


    
  def test_add_31_21(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4

    # Test most common case where type(lcm(D1, D2)) = 51 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1], [0, 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[1, 0, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + 1, z4^3 + z4 + 1, z4^3 + z4 + 1, 1], [z4^3 + z4^2, z4^3 + z4 + 1, z4^2, 0, 1], [z4 + 1, z4^3 + z4^2 + 1, z4^3 + z4^2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3, z4^3 + z4^2, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4, z4^2 + 1, 1], [z4^3 + z4^2 + 1, 0, z4^3 + z4 + 1, 0, 1], [z4^3 + z4^2 + 1, z4^2 + z4, z4^2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 2, 0, 0, 1], [1, 2, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1, 1], [0, 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[1, 0, 1, 1], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_3_3, [[z3^2 + z3, z3^2 + 1, z3^2 + z3 + 1, 1], [2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 1, 0, 1], [z3^2, z3^2 + z3, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2 + z3 + 2, z3^2 + z3, 1], [z3^2 + 2*z3, 2*z3^2 + 2*z3 + 2, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[z3^2, 0, z3^2 + z3 + 2, 1], [2*z3^2 + 2*z3 + 1, z3 + 1, z3^2 + 2, 0, 1], [2*z3^2 + 2, z3^2 + z3 + 1, 2*z3^2 + 2*z3, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31, [[6, 20, 7, 1], [6, 19, 16, 0, 1], [4, 16, 24, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[23, 26, 1], [8, 4, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[16, 0, 20, 1], [11, 18, 23, 0, 1], [10, 0, 22, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_31_2, [[3*z2 + 19, 9*z2 + 22, 6*z2 + 9, 1], [10*z2 + 18, 10*z2 + 6, 24*z2 + 3, 0, 1], [22*z2, 4*z2 + 21, 22*z2 + 16, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[7*z2 + 10, 13*z2 + 19, 1], [2*z2 + 9, 28*z2 + 10, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[25*z2 + 14, 11*z2 + 22, 10*z2 + 21, 1], [z2 + 15, 13*z2 + 27, 8*z2 + 28, 0, 1], [20*z2 + 27, 26*z2, 11*z2 + 20, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_1009, [[505, 854, 752, 1], [933, 33, 629, 0, 1], [167, 606, 248, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[847, 396, 1], [816, 514, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[32, 831, 986, 1], [91, 1000, 91, 0, 1], [703, 566, 373, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 52
    D1 = C34CurveDivisor(C_31, [[13, 16, 1, 1], [28, 26, 0, 0, 1], [14, 3, 8, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[10, 1, 1], [13, 5, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[1, 1], [29, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 53
    D1 = C34CurveDivisor(C_1009, [[977, 497, 31, 1], [586, 808, 458, 0, 1], [642, 722, 719, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[458, 559, 1], [918, 432, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[284, 92, 1], [108, 541, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 54
    D1 = C34CurveDivisor(C_31, [[19, 13, 1, 1], [2, 14, 10, 0, 1], [12, 14, 3, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[5, 24, 1], [14, 20, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[10, 1], [24, 0, 21, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 41
    D1 = C34CurveDivisor(C_1009, [[374, 880, 892, 1], [719, 789, 169, 0, 1], [476, 712, 226, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[199, 994, 1], [471, 36, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[228, 949, 149, 1], [875, 287, 57, 0, 1], [273, 517, 787, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test case where type(lcm(D1, D2)) = 42
    D1 = C34CurveDivisor(C_1009, [[392, 897, 0, 1], [486, 893, 48, 0, 1], [675, 299, 814, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[670, 714, 1], [392, 897, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[355, 243, 1], [375, 689, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 43
    C = C34Curve(GF(37), [4, 23, 16, 2, 30, 0, 0, 1, 0])
    D1 = C34CurveDivisor(C, [[5, 11, 0, 1], [13, 15, 28, 0, 1], [6, 9, 33, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[24, 6, 1], [5, 11, 0, 1], []])
    D3 = C34CurveDivisor(C, [[26, 32, 1], [30, 3, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    C = C34Curve(GF(67), [46, 58, 50, 64, 63, 5, 66, 51, 24])
    D1 = C34CurveDivisor(C, [[59, 65, 0, 1], [19, 12, 63, 0, 1], [5, 65, 63, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[43, 49, 1], [59, 65, 0, 1], []])
    D3 = C34CurveDivisor(C, [[52, 20, 1], [4, 4, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where D2 < D1
    D1 = C34CurveDivisor(C_1009, [[155, 847, 246, 1], [25, 479, 425, 0, 1], [706, 729, 91, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[50, 299, 1], [972, 950, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[275, 187, 7, 1], [346, 592, 934, 0, 1], [361, 959, 988, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)



  def test_add_31_22(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 51 under various fields
    D1 = C34CurveDivisor(C_2, [[1, 0, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[1, 1], [0, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3, 1, z4^3 + z4^2 + z4, 1], [1, z4^3 + 1, 1, 0, 1], [z4^3 + z4^2, z4^3 + z4 + 1, z4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, 1], [z4^3 + z4^2 + z4, 0, 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_2_4, [[0, z4^3 + z4^2 + z4, z4^2 + z4, 1], [z4^3 + z4^2 + z4, z4^3 + z4^2 + 1, z4^2, 0, 1], [1, z4^3 + z4 + 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[0, 2, 1, 1], [1, 0, 0, 0, 1], [1, 2, 0, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[2, 1], [0, 0, 2, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 2, 0, 1], [0, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 0, 0, 1], [z3^2 + 2*z3 + 1, z3 + 2, 2*z3^2, 0, 1], [2*z3^2 + z3, 2*z3^2 + 1, z3^2 + z3 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3, 1], [z3^2 + 2, 0, z3^2 + 2*z3 + 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + 2*z3 + 2, 2*z3^2 + 2*z3 + 2, 2*z3^2 + z3, 1], [z3^2 + 2, 0, z3^2 + z3 + 1, 0, 1], [2*z3, 2*z3^2 + 2, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[13, 20, 30, 1], [29, 9, 0, 0, 1], [29, 29, 27, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[11, 1], [26, 0, 27, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 11, 12, 1], [18, 11, 4, 0, 1], [3, 28, 16, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[11*z2 + 28, 20*z2 + 27, 18*z2 + 27, 1], [18*z2 + 21, 3*z2 + 28, 27*z2 + 23, 0, 1], [z2 + 16, 28*z2 + 23, 9*z2 + 4, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 28, 1], [9*z2 + 12, 0, 7*z2 + 1, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31_2, [[10*z2 + 22, 29*z2 + 13, 7*z2 + 6, 1], [11, 22*z2 + 20, 25*z2 + 9, 0, 1], [19*z2 + 7, 17*z2 + 9, 16*z2 + 26, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_1009, [[35, 840, 550, 1], [118, 47, 846, 0, 1], [340, 960, 370, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[436, 1], [327, 0, 276, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[529, 488, 1, 1], [538, 172, 205, 0, 1], [283, 416, 199, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # TODO: Test case where type(lcm(D1, D2)) = 52 (This case is impossible?)
    # Test case where type(lcm(D1, D2)) = 53
    D1 = C34CurveDivisor(C_1009, [[830, 401, 629, 1], [626, 655, 620, 0, 1], [747, 11, 800, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[652, 1], [460, 0, 399, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[708, 446, 1], [467, 369, 0, 1], []])
    self.assertEqual(D1 + D2, D3)
    
    # Test case where type(lcm(D1, D2)) = 54
    D1 = C34CurveDivisor(C_31, [[30, 20, 0, 1], [24, 20, 26, 0, 1], [10, 25, 28, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[25, 1], [28, 0, 14, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[26, 1], [11, 0, 28, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 41
    D1 = C34CurveDivisor(C_1009, [[706, 1006, 543, 1], [89, 400, 110, 0, 1], [245, 315, 86, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[328, 1], [555, 0, 277, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[962, 541, 818, 1], [114, 164, 50, 0, 1], [307, 375, 409, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 42
    D1 = C34CurveDivisor(C_31, [[16, 17, 0, 1], [2, 2, 1, 0, 1], [14, 27, 10, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[1, 1], [10, 0, 24, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[9, 14, 1], [16, 17, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 43
    D1 = C34CurveDivisor(C_31, [[2, 4, 0, 1], [24, 27, 25, 0, 1], [6, 14, 14, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[10, 1], [24, 0, 21, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_31, [[10, 13, 5, 1], [14, 20, 10, 0, 1], [13, 0, 16, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where D2 < D1
    D1 = C34CurveDivisor(C_1009, [[183, 407, 0, 1], [632, 238, 469, 0, 1], [234, 33, 78, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[469, 1], [901, 0, 78, 0, 0, 1], []])
    D3 = C34CurveDivisor(C_1009, [[199, 733, 687, 1], [801, 37, 659, 0, 1], [954, 472, 755, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)



  def test_add_31_31(self) :
    C_2, C_3, C_31, C_1009 = self.C_2, self.C_3, self.C_31, self.C_1009
    C_2_4, C_3_3, C_31_2 = self.C_2_4, self.C_3_3, self.C_31_2
    z2, z3, z4 = self.z2, self.z3, self.z4
    
    # Test most common case where type(lcm(D1, D2)) = 61 under various fields
    D1 = C34CurveDivisor(C_2, [[0, 1, 0, 1], [0, 1, 0, 0, 1], [1, 1, 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2, [[0, 1, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1]])
    D3 = C34CurveDivisor(C_2, [[1, 0, 0, 1], [0, 0, 1, 0, 1], [1, 1, 1, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    D1 = C34CurveDivisor(C_2_4, [[z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + z4 + 1, 1, 1], [z4^3, z4^2 + z4 + 1, z4^3 + z4^2, 0, 1], [0, z4^3 + z4, z4^3 + z4^2 + z4 + 1, 0, 0, 1]])
    D2 = C34CurveDivisor(C_2_4, [[z4^3 + z4 + 1, z4^3 + z4, z4^2 + z4 + 1, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + 1, z4^2 + z4 + 1, 0, 1], [z4^3 + z4^2 + z4 + 1, z4^3 + z4^2 + 1, z4^3 + z4 + 1, 0, 0, 1]])
    D3 = C34CurveDivisor(C_2_4, [[z4^2 + z4, z4^3 + z4^2 + z4 + 1, z4 + 1, 1], [z4^2 + 1, z4 + 1, z4, 0, 1], [z4, 1, z4^2 + z4, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3, [[0, 1, 0, 1], [0, 2, 0, 0, 1], [1, 1, 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3, [[1, 1, 2, 1], [2, 1, 0, 0, 1], [1, 2, 0, 0, 0, 1]])
    D3 = C34CurveDivisor(C_3, [[0, 2, 0, 1], [0, 0, 2, 0, 1], [0, 0, 2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_3_3, [[2*z3^2, 2*z3^2 + 1, z3^2, 1], [2*z3, 2*z3^2 + 2*z3 + 2, 2*z3^2, 0, 1], [z3 + 2, z3 + 2, 2*z3^2 + 2*z3 + 2, 0, 0, 1]])
    D2 = C34CurveDivisor(C_3_3, [[z3^2, z3^2 + 2, 0, 1], [z3^2 + 2, z3^2 + 2*z3 + 2, 2*z3^2 + z3 + 2, 0, 1], [z3^2 + 2, z3 + 1, z3^2, 0, 0, 1]])
    D3 = C34CurveDivisor(C_3_3, [[2*z3^2 + z3, 2*z3^2 + 2*z3 + 2, z3^2 + 1, 1], [2*z3^2 + 2*z3, 2*z3^2 + 2, z3^2 + 2*z3 + 1, 0, 1], [2*z3^2, z3^2 + z3 + 2, z3^2, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31, [[26, 4, 2, 1], [17, 24, 22, 0, 1], [0, 23, 18, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[9, 21, 4, 1], [0, 2, 4, 0, 1], [20, 2, 26, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31, [[28, 9, 21, 1], [15, 20, 21, 0, 1], [19, 6, 23, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_31_2, [[5*z2 + 18, 10*z2 + 19, 9*z2 + 13, 1], [19*z2 + 30, 11*z2 + 3, 12*z2 + 16, 0, 1], [6*z2 + 24, 10*z2 + 6, 10*z2 + 5, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31_2, [[24*z2 + 26, z2 + 29, 30*z2 + 6, 1], [21*z2 + 2, 14*z2 + 28, 8, 0, 1], [5*z2 + 6, 25*z2 + 9, 10*z2 + 27, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31_2, [[7*z2 + 1, 5, 11*z2 + 3, 1], [16*z2 + 26, 14*z2 + 9, 12*z2 + 26, 0, 1], [16*z2 + 9, 13*z2 + 3, z2 + 18, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    D1 = C34CurveDivisor(C_1009, [[351, 509, 773, 1], [395, 206, 757, 0, 1], [450, 72, 961, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[1007, 845, 679, 1], [452, 771, 469, 0, 1], [165, 413, 223, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[877, 43, 18, 1], [191, 561, 937, 0, 1], [1000, 380, 341, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 62
    D1 = C34CurveDivisor(C_31, [[4, 13, 5, 1], [6, 26, 22, 0, 1], [13, 14, 23, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[23, 8, 25, 1], [18, 21, 15, 0, 1], [7, 28, 10, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31, [[10, 1], [25, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 63
    D1 = C34CurveDivisor(C_1009, [[386, 75, 575, 1], [207, 322, 985, 0, 1], [907, 188, 32, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[481, 885, 975, 1], [914, 637, 160, 0, 1], [299, 353, 147, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[673, 724, 1], [829, 971, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 64
    D1 = C34CurveDivisor(C_1009, [[135, -399, -106, 1], [298, -265, -212, 0, 1], [-331, -432, 498, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[234, -479, 421, 1], [-242, -379, -243, 0, 1], [-284, 459, 305, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[342, 1], [-500, 0, -370, 0, 0, 1], []])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 65
    D1 = C34CurveDivisor(C_1009, [[361, -264, -47, 1], [-159, 214, -297, 0, 1], [370, -390, 191, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[361, -264, -47, 1], [354, 425, -344, 0, 1], [494, -467, 25, 0, 0, 1]])
    D3 = C_1009.zero_divisor()
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 51
    D1 = C34CurveDivisor(C_1009, [[0, 883, 370, 1], [545, 467, 772, 0, 1], [341, 43, 639, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[669, 291, 248, 1], [966, 973, 427, 0, 1], [514, 105, 506, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[651, 2, 608, 1], [915, 260, 108, 0, 1], [262, 979, 567, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    # Test case where type(lcm(D1, D2)) = 52
    C = C34Curve(GF(997), [44, 762, 377, 376, 677, 281, 243, 890, 703])
    D1 = C34CurveDivisor(C, [[789, 913, 741, 1], [817, 236, 194, 0, 1], [543, 428, 696, 0, 0, 1]])
    D2 = C34CurveDivisor(C, [[621, 466, 432, 1], [429, 984, 691, 0, 1], [454, 209, 16, 0, 0, 1]])
    D3 = C34CurveDivisor(C, [[839, 891, 1], [376, 526, 0, 1], []])    
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 53
    D1 = C34CurveDivisor(C_31, [[3, 3, 3, 1], [25, 5, 1, 0, 1], [9, 14, 26, 0, 0, 1]])
    D2 = C34CurveDivisor(C_31, [[14, 11, 25, 1], [25, 5, 1, 0, 1], [30, 24, 25, 0, 0, 1]])
    D3 = C34CurveDivisor(C_31, [[19, 10, 2, 1], [22, 11, 4, 0, 1], [23, 11, 24, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 54
    D1 = C34CurveDivisor(C_1009, [[372, 989, 697, 1], [848, 907, 18, 0, 1], [691, 332, 826, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[372, 989, 697, 1], [40, 105, 471, 0, 1], [187, 925, 944, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[107, 209, 0, 1], [644, 102, 402, 0, 1], [372, 744, 267, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)
    
    # Test case where type(lcm(D1, D2)) = 41
    D1 = C34CurveDivisor(C_1009, [[146, 655, 26, 1], [25, 32, 557, 0, 1], [784, 801, 189, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[545, 621, 740, 1], [148, 985, 724, 0, 1], [172, 489, 687, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[487, 436, 734, 1], [40, 447, 302, 0, 1], [901, 632, 9, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 42
    D1 = C34CurveDivisor(C_1009, [[766, 943, 0, 1], [647, 996, 804, 0, 1], [364, 718, 185, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[766, 943, 0, 1], [647, 996, 804, 0, 1], [463, 25, 762, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[294, 211, 592, 1], [989, 620, 871, 0, 1], [981, 435, 61, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

    # Test case where type(lcm(D1, D2)) = 43
    D1 = C34CurveDivisor(C_1009, [[38, 620, 623, 1], [960, 86, 938, 0, 1], [424, 447, 594, 0, 0, 1]])
    D2 = C34CurveDivisor(C_1009, [[38, 620, 623, 1], [967, 296, 350, 0, 1], [412, 87, 593, 0, 0, 1]])
    D3 = C34CurveDivisor(C_1009, [[33, 779, 445, 1], [607, 232, 156, 0, 1], [970, 182, 1000, 0, 0, 1]])
    self.assertEqual(D1 + D2, D3)

