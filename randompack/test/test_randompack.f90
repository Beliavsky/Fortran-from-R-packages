program test_randompack
   use iso_fortran_env, only : int8, int64
   use ieee_arithmetic, only : ieee_is_finite
   use randompack, only : dp, randompack_rng, randompack_rng_type, randompack_snapshot
   use randompack_engines, only : engine_state_type, engine_name, seed_engine, next_u64
   use randompack_openlibm, only : openlibm_log, openlibm_log1p, openlibm_exp
   implicit none
   integer(int64), parameter :: expected(4, 14) = reshape([ &
      int(z'51008A03CF775CA2', int64), int(z'D85811148A5E1E45', int64), &
      int(z'2A4E4AF819C90612', int64), int(z'7119503073C28FBD', int64), &
      int(z'FEBF80CA4B5DE5A7', int64), int(z'CF9F5818D3874B00', int64), &
      int(z'5B1203D276B6C6B3', int64), int(z'ABF1F42CC22C6E33', int64), &
      int(z'39067EC9BF7F330B', int64), int(z'59067EC9BF7F330B', int64), &
      int(z'79067EC9BF7F330B', int64), int(z'99067EC9BF7F330B', int64), &
      int(z'51008A03CF775CA2', int64), int(z'7903BA068B6C61F6', int64), &
      int(z'39C2AA07B61A524A', int64), int(z'3D894C0FAC39304B', int64), &
      int(z'FEBF80CA4B5DE5A7', int64), int(z'A7AF855C48087BAF', int64), &
      int(z'0AE0C40DA269E193', int64), int(z'EAB0FE1D3B66A18F', int64), &
      int(z'69B673137C0143EF', int64), int(z'58A4D47D8CF84C5C', int64), &
      int(z'CD51C2B57C030DEA', int64), int(z'F132370EC5BE0C89', int64), &
      int(z'04480F4CC57765BD', int64), int(z'25DF737F0D15FA41', int64), &
      int(z'CA66115DEBEDBFBD', int64), int(z'FB0365BF16C0F314', int64), &
      int(z'28405A10E39BE950', int64), int(z'78E185654D7C69AB', int64), &
      int(z'0707393CEB3E13AC', int64), int(z'27C9C886EC70B576', int64), &
      int(z'39067EC9BF7F330B', int64), int(z'1EF8AA439CD7FA00', int64), &
      int(z'C39B0114BB92BAA8', int64), int(z'ACA2ABD3A56254DE', int64), &
      int(z'3B45FB31E53E3575', int64), int(z'E1CFF273D1604B4C', int64), &
      int(z'50BD1285FAD1183D', int64), int(z'AD40BAD2A0004912', int64), &
      int(z'B317FF3D3B1A6D71', int64), int(z'AAACA15ADE038ED9', int64), &
      int(z'2964538838F252A5', int64), int(z'4CC98A582B4C50B9', int64), &
      int(z'A4D57C2FA8B1CC5B', int64), int(z'0749A591B47BE11A', int64), &
      int(z'389B04BCE0D736D1', int64), int(z'DE5F0D9EA283BDEF', int64), &
      int(z'2A41D50E8D062A1C', int64), int(z'4F45ADC57F4C6438', int64), &
      int(z'C72E9F7A557BD130', int64), int(z'80562F86ABF854A0', int64), &
      int(z'0EEE681B767BA9E5', int64), int(z'933296D72A48713B', int64), &
      int(z'24A4C8511E557288', int64), int(z'AE9D6C7950B1315D', int64)], [4, 14])
   integer, parameter :: expected_int(64) = [ &
      15, 15, 3, -2, 4, 1, 6, 5, -4, 16, 0, -5, 11, 16, 13, -2, &
      6, 10, -3, 16, 1, 16, 16, 11, 12, 16, 7, 16, 10, 13, 3, 0, &
      11, -2, 10, 10, 5, 7, 8, 15, 1, 9, 13, 15, -5, -2, 1, -5, &
      5, 1, 11, -2, 6, 1, -3, -2, 14, -3, -3, 2, 9, -5, -4, 3]
   integer, parameter :: expected_perm(25) = [ &
      23, 25, 3, 9, 22, 10, 20, 1, 24, 4, 18, 17, 11, &
      6, 15, 8, 19, 12, 21, 7, 5, 14, 2, 13, 16]
   integer, parameter :: expected_sample(10) = [7, 8, 32, 20, 6, 25, 10, 30, 36, 40]
   integer(int64), parameter :: expected_norm_bits(10) = [ &
      int(z'BFD3FCC3FB80770E', int64), int(z'BFF6F7408FBED88F', int64), &
      int(z'BFDA0F916BCA5109', int64), int(z'BFD6DDA97212829E', int64), &
      int(z'3FED316776830DD0', int64), int(z'BFF1C875313E7EB7', int64), &
      int(z'BFE09B0C128356B6', int64), int(z'3FE9439421C89F54', int64), &
      int(z'BFEDECF76F192B5B', int64), int(z'BFDC08F62EF8BB67', int64)]
   integer(int64), parameter :: expected_exp_bits(10) = [ &
      int(z'3FD2805B07C7E90C', int64), int(z'3FFC7B52F9D2358A', int64), &
      int(z'3FE035E2B3C83FAA', int64), int(z'3FD81D2A8F675BA0', int64), &
      int(z'3FEE342E47C61C41', int64), int(z'3FF68B92A3F2E72C', int64), &
      int(z'3FD9330B9FDC0FD7', int64), int(z'3FE4C6E929ECA278', int64), &
      int(z'3FF36739BC1BDEAB', int64), int(z'3FEF2C6822060741', int64)]
   integer(int64), parameter :: expected_gamma_hi_bits(17) = [ &
      int(z'4000B1106C8BFE2A', int64), int(z'3FE99529BDF8A3C2', int64), &
      int(z'3FFF1ED4D2AB7BEC', int64), int(z'400026188A0C6F28', int64), &
      int(z'401244D9A56CDD01', int64), int(z'3FF16E780BF3E268', int64), &
      int(z'3FFC9717475466A9', int64), int(z'401108ECE1C9C524', int64), &
      int(z'3FF45F56F667E36E', int64), int(z'3FFE6854277E348F', int64), &
      int(z'400EBD23A9874E08', int64), int(z'4001A075D295953A', int64), &
      int(z'40034A08051AA5A7', int64), int(z'40104FBA5EA3CF66', int64), &
      int(z'40123C0D7446BF53', int64), int(z'40218EA5ACDFDF9B', int64), &
      int(z'400D159762D46D79', int64)]
   integer(int64), parameter :: expected_gamma_lo_bits(17) = [ &
      int(z'40000EA936E492AE', int64), int(z'3FC6B7784EAA037B', int64), &
      int(z'3FD8733C0DC81300', int64), int(z'3FFAB26D0E912E95', int64), &
      int(z'4008BBC281887D3D', int64), int(z'3FD774AF3B8F8C4A', int64), &
      int(z'3FEA4BF9DCA0B6C9', int64), int(z'400B92C246F3C5BC', int64), &
      int(z'3FD658560B3A1062', int64), int(z'3FE00D370C937EB5', int64), &
      int(z'3FFBC887738C6BC9', int64), int(z'3F5B699FC9BB63F8', int64), &
      int(z'4002114231849B0F', int64), int(z'400F9FC4CA02AA11', int64), &
      int(z'400CBCE43A05417C', int64), int(z'3FF92CA1D05A65CC', int64), &
      int(z'4000C6EA67BA3651', int64)]
   integer(int64), parameter :: expected_beta_bits(17) = [ &
      int(z'3FC08AB35291E2E7', int64), int(z'3FAE32D8714CB77D', int64), &
      int(z'3FC3132D8D49B22E', int64), int(z'3FCC7828DBDCB4AF', int64), &
      int(z'3FD6843EE430D368', int64), int(z'3FB10F56F3839BDD', int64), &
      int(z'3FD1C2C5F17C61A9', int64), int(z'3FDAB42F2D020942', int64), &
      int(z'3FB5AA1876B5B280', int64), int(z'3FCA0185CC51358F', int64), &
      int(z'3FD58CE154B5C480', int64), int(z'3FCE470C30652872', int64), &
      int(z'3FD8A6477EF10691', int64), int(z'3FD7A046A769A867', int64), &
      int(z'3FE6A7EC937CE454', int64), int(z'3FDAD5C483C31851', int64), &
      int(z'3FE2343693CF77B1', int64)]
   integer(int64), parameter :: expected_t_bits(17) = [ &
      int(z'BFD3F0EF04F57D96', int64), int(z'BFFC7C00913D0F7F', int64), &
      int(z'BFD27A0C174C2C16', int64), int(z'BFE33A5EFF3986A1', int64), &
      int(z'3FE1BB4BE2DA9B3A', int64), int(z'BFF0980BD7BC3664', int64), &
      int(z'BFE2BC1C27944729', int64), int(z'40018D2EF756C063', int64), &
      int(z'BFF231911284F934', int64), int(z'BFDD62CD5FE76889', int64), &
      int(z'3FE20B226F7E184B', int64), int(z'BFCF9226A672A124', int64), &
      int(z'BFBB710290A1B6A6', int64), int(z'3FF390CCC177CD13', int64), &
      int(z'3FF4965CEF781176', int64), int(z'40042BCDDCFB5B0A', int64), &
      int(z'3FD86C107114EDCB', int64)]
   integer(int64), parameter :: expected_f_bits(17) = [ &
      int(z'3FDA1007495085B1', int64), int(z'3FC9FC02349F7EAB', int64), &
      int(z'3FDEC3FB4815F5DF', int64), int(z'3FE875436AC96E10', int64), &
      int(z'3FF51C67313E35C3', int64), int(z'3FCC086EAED97639', int64), &
      int(z'3FF0597656D29E7B', int64), int(z'3FFB968B2F2F276C', int64), &
      int(z'3FD1908F07619E48', int64), int(z'3FE60F0407CD8A63', int64), &
      int(z'3FF407784BDBE7E1', int64), int(z'3FEA2BA3FB759F61', int64), &
      int(z'3FF940A33612FDD0', int64), int(z'3FF6D6DDB7A9EB2D', int64), &
      int(z'4015940511906F7F', int64), int(z'3FFACA8001DE55AC', int64), &
      int(z'4008BBBC881AC661', int64)]
   integer(int64), parameter :: expected_lognormal_bits(17) = [ &
      int(z'3FEBE5953A1D2D09', int64), int(z'3FC72983CF5C5648', int64), &
      int(z'3FE86D1A51E5A56E', int64), int(z'3FEA31C333A5219E', int64), &
      int(z'40135D8A3362CEF7', int64), int(z'3FD239F0BF8E6A7B', int64), &
      int(z'3FE4E3AB3EC9C4CC', int64), int(z'40104E9676B8E5EA', int64), &
      int(z'3FD753DB8CC72AE4', int64), int(z'3FE764CEB4F5CE31', int64), &
      int(z'40098278CA6AD094', int64), int(z'3FEEFEF5224E63CC', int64), &
      int(z'3FF2866EBFD1D2D3', int64), int(z'400D5F8C613E4DA6', int64), &
      int(z'401346643CEFA56D', int64), int(z'403DBCA22FF808CB', int64), &
      int(z'40066DF81D2E807A', int64)]
   integer(int64), parameter :: expected_gumbel_bits(17) = [ &
      int(z'BFE33145CCD5572E', int64), int(z'3FE258D3D8EE68D5', int64), &
      int(z'BFF76A7A51066352', int64), int(z'BFE373AFBACB1BE6', int64), &
      int(z'4010CA65F4C00983', int64), int(z'3FF87A699C9288FB', int64), &
      int(z'401002138BB67EE9', int64), int(z'BFD511F6090F129B', int64), &
      int(z'3FF761D52C50B59F', int64), int(z'4007DB190A92EAC0', int64), &
      int(z'400A5C07ACA5B314', int64), int(z'BFF96E284C115A09', int64), &
      int(z'BFE572E9D414DC06', int64), int(z'BFE1213561EE3C86', int64), &
      int(z'3FB7FC6EFE6419FF', int64), int(z'3FC81467BBCDFA59', int64), &
      int(z'4000B00605EF20EF', int64)]
   integer(int64), parameter :: expected_pareto_bits(17) = [ &
      int(z'3FFED7C6A2B8EF54', int64), int(z'400D7D4795724B81', int64), &
      int(z'4000F36E48D80C42', int64), int(z'400005527F0DC258', int64), &
      int(z'4004801EF6F0F8FE', int64), int(z'400918980E9E0C4B', int64), &
      int(z'400023AD731EE445', int64), int(z'40020932C301F6BB', int64), &
      int(z'40070AD3AF399593', int64), int(z'4004C5B81247AA0B', int64), &
      int(z'40038DD7DD9B3E36', int64), int(z'3FFE7589D715F0F1', int64), &
      int(z'3FFC25078C6BFAE6', int64), int(z'40042ABFEAE8005C', int64), &
      int(z'4003EED9F877914E', int64), int(z'401D17BC8EF236E7', int64), &
      int(z'4001843CB9F8DC05', int64)]
   integer(int64), parameter :: expected_weibull_bits(17) = [ &
      int(z'3FF03118F4EFDAC0', int64), int(z'400795B037B1F812', int64), &
      int(z'3FF685900D1ACD62', int64), int(z'3FF2EC2B326F361B', int64), &
      int(z'40003D1C720546D3', int64), int(z'40048E1503583B44', int64), &
      int(z'3FF36B45C43BAD36', int64), int(z'3FFA0FD6D91B02B0', int64), &
      int(z'4002D170800D9CE6', int64), int(z'40008B1A6E9D1AE5', int64), &
      int(z'3FFE38CFC96C3751', int64), int(z'3FEE71FB17F9BC03', int64), &
      int(z'3FDE16914D76EBDD', int64), int(z'3FFFB5200E16D4E0', int64), &
      int(z'3FFF26D19EDD4B38', int64), int(z'401115B602A091AF', int64), &
      int(z'3FF87081C4B8CFFD', int64)]
   integer(int64), parameter :: expected_skew_normal_bits(17) = [ &
      int(z'3FA95CE141DD2EA8', int64), int(z'BFEB882894279443', int64), &
      int(z'BFF916A86277BFF0', int64), int(z'BFF5BD50987CFA26', int64), &
      int(z'BFFD1ED40CF41F8C', int64), int(z'BFE2CCDDCB6F25D6', int64), &
      int(z'BFBACFA74F043A68', int64), int(z'C000ADC5BFF028E0', int64), &
      int(z'BFE220C0EC5436F2', int64), int(z'3FC8B29E116217E5', int64), &
      int(z'3FBCC6B2139ABFE4', int64), int(z'3FD031C99532CAFD', int64), &
      int(z'3FAD3B0B418FE5D8', int64), int(z'BFEDC5A01117318F', int64), &
      int(z'BFD0F6F313FDAD58', int64), int(z'3FF111E3667D5CF6', int64), &
      int(z'BFEFF44D7529136B', int64)]
   integer(int64), parameter :: expected_mvn_full_bits(15) = [ &
      int(z'3FE6019E023FC479', int64), int(z'BFDBDD023EFB623C', int64), &
      int(z'3FE2F8374A1AD77C', int64), int(z'3FE4912B46F6BEB1', int64), &
      int(z'3FFE98B3BB4186E8', int64), int(z'C005004683D8ECF5', int64), &
      int(z'C000A63555031E92', int64), int(z'BF6775214B633900', int64), &
      int(z'C00320D4994F1EB7', int64), int(z'BFF61B4D2B1F02EC', int64), &
      int(z'3FEFA4423C1F1EB4', int64), int(z'BFA02EFD631190E0', int64), &
      int(z'3FDEED152985A714', int64), int(z'3FF23F5A06729958', int64), &
      int(z'3FF997C82ABDBBF4', int64)]
   integer(int64), parameter :: expected_mvn_singular_bits(15) = [ &
      int(z'3FF1ABC52EDB8453', int64), int(z'3FB27D835D99A928', int64), &
      int(z'BFE46C960ED5780F', int64), int(z'3FEDA74218BC48F5', int64), &
      int(z'3FF0CBEBDFCF592D', int64), int(z'C0006C1087EFF47F', int64), &
      int(z'C00C412045E4822A', int64), int(z'BFE0B6742412A845', int64), &
      int(z'BFFF9BECF2758EEE', int64), int(z'3FE7AF326AA54A64', int64), &
      int(z'3FD34E907BEE6D5C', int64), int(z'C001AD342AF7B4E1', int64), &
      int(z'3FB6E7AE68BEFD40', int64), int(z'3FC9BDA0CF44AC74', int64), &
      int(z'400851C28A90FF30', int64)]
   integer(int64), parameter :: expected_mvn_full_post_bits(6) = [ &
      int(z'3FD75445010199E8', int64), int(z'3FE9503FC56A34D8', int64), &
      int(z'3FC40CEEF455C250', int64), int(z'3FED1999A8B085CC', int64), &
      int(z'3FB0A3F83E581290', int64), int(z'3FD4CB361364A82C', int64)]
   integer(int64), parameter :: expected_mvn_singular_post_bits(6) = [ &
      int(z'3FED2C59B154C5E8', int64), int(z'3F93811220A94C00', int64), &
      int(z'3FC2275E08A4A3B8', int64), int(z'3FC6028B9D23C5C0', int64), &
      int(z'3FD59F0694C6EB34', int64), int(z'3FD75445010199E8', int64)]
   integer(int64), parameter :: expected_openlibm_bits(6) = [ &
      int(z'BFE62E42FEFA39EF', int64), int(z'C03BA18C0CB51D2A', int64), &
      int(z'BD719799812DF3BD', int64), int(z'3FE368B2FC6F960A', int64), &
      int(z'0000000000000001', int64), int(z'7F0D945DF4F8EC8E', int64)]
   integer(int8), parameter :: expected_raw(31) = [ &
      int(z'1B', int8), int(z'5E', int8), int(z'88', int8), int(z'A3', int8), &
      int(z'E7', int8), int(z'61', int8), int(z'21', int8), int(z'11', int8), &
      int(z'EE', int8), int(z'E9', int8), int(z'74', int8), int(z'98', int8), &
      int(z'36', int8), int(z'D2', int8), int(z'ED', int8), int(z'C5', int8), &
      int(z'3B', int8), int(z'80', int8), int(z'B6', int8), int(z'0A', int8), &
      int(z'62', int8), int(z'28', int8), int(z'84', int8), int(z'E3', int8), &
      int(z'52', int8), int(z'E8', int8), int(z'16', int8), int(z'6A', int8), &
      int(z'AC', int8), int(z'8E', int8), int(z'D5', int8)]
   type(engine_state_type) :: state
   type(randompack_rng_type) :: a, b
   type(randompack_snapshot) :: snap
   real(dp) :: x(32), y(32), mv(16, 2), sigma(2, 2), mu(2)
   real(dp) :: mv3(5, 3), sigma3(3, 3), mu3(3), post6(6)
   integer(int64) :: word, delta(2)
   integer(int8) :: bytes(31)
   integer :: e, i, info, ints(64), j, k, p(25), s(10)

   do e = 1, 14
      state = engine_state_type()
      state%engine = e
      call seed_engine(state, 123, info=info)
      call require(info == 0, 'seed '//trim(engine_name(e)))
      do i = 1, 4
         call next_u64(state, word)
         call require(word == expected(i, e), 'engine vector '//trim(engine_name(e)))
      end do
   end do

   a = randompack_rng('pcg64', seed=321)
   b = a%duplicate()
   call a%unif(x, info=info)
   call b%unif(y, info=info)
   call require(all(x == y), 'duplicate stream')
   snap = a%serialize()
   call a%unif(x, info=info)
   call b%deserialize(snap)
   call b%unif(y, info=info)
   call require(all(x == y), 'snapshot restore')

   a = randompack_rng('pcg64', seed=123)
   b = randompack_rng('pcg64', seed=123)
   call a%jump(80, info)
   delta = [0_int64, 65536_int64]
   call b%advance(delta, info)
   call a%unif(x(1:4), info=info)
   call b%unif(y(1:4), info=info)
   call require(all(x(1:4) == y(1:4)), 'pcg jump equals advance')

   call require(transfer(openlibm_log(0.5_dp), 0_int64) == expected_openlibm_bits(1), 'openlibm log reference')
   call require(transfer(openlibm_log1p(-0.999999999999_dp), 0_int64) == expected_openlibm_bits(2), &
      'openlibm log1p near -1 reference')
   call require(transfer(openlibm_log1p(-1.0e-12_dp), 0_int64) == expected_openlibm_bits(3), &
      'openlibm log1p small reference')
   call require(transfer(openlibm_exp(-0.5_dp), 0_int64) == expected_openlibm_bits(4), 'openlibm exp reference')
   call require(transfer(openlibm_exp(-745.0_dp), 0_int64) == expected_openlibm_bits(5), 'openlibm exp subnormal reference')
   call require(transfer(openlibm_exp(700.0_dp), 0_int64) == expected_openlibm_bits(6), 'openlibm exp large reference')

   a = randompack_rng('pcg64', seed=9071, bitexact=.true.)
   call a%normal(x(1:1), info=info)
   call require(info == 0 .and. transfer(x(1), 0_int64) == int(z'400EF9C97BE4D17D', int64), &
      'upstream normal Ziggurat tail bit vector')
   a = randompack_rng('pcg64', seed=753, bitexact=.true.)
   call a%exp(x(1:1), info=info)
   call require(info == 0 .and. transfer(x(1), 0_int64) == int(z'4020DC22F1CF89A2', int64), &
      'upstream exponential Ziggurat tail bit vector')

   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%normal(x(1:10), info=info)
   call require(info == 0, 'upstream normal status')
   do i = 1, 10
      call require(transfer(x(i), 0_int64) == expected_norm_bits(i), 'upstream normal bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%exp(x(1:10), info=info)
   call require(info == 0, 'upstream exponential status')
   do i = 1, 10
      call require(transfer(x(i), 0_int64) == expected_exp_bits(i), 'upstream exponential bit vector')
   end do

   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%gamma(x(1:17), 2.5_dp, 1.2_dp, info)
   call require(info == 0, 'upstream gamma shape>=1 status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_gamma_hi_bits(i), &
         'upstream gamma shape>=1 bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%gamma(x(1:17), 0.7_dp, 2.0_dp, info)
   call require(info == 0, 'upstream gamma shape<1 status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_gamma_lo_bits(i), &
         'upstream gamma shape<1 bit vector')
   end do

   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%beta(x(1:17), 2.0_dp, 5.0_dp, info)
   call require(info == 0, 'upstream beta status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_beta_bits(i), 'upstream beta bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%t(x(1:17), 7.0_dp, info)
   call require(info == 0, 'upstream Student-t status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_t_bits(i), 'upstream Student-t bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%f(x(1:17), 5.0_dp, 11.0_dp, info)
   call require(info == 0, 'upstream F status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_f_bits(i), 'upstream F bit vector')
   end do

   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%lognormal(x(1:17), 0.3_dp, 1.4_dp, info)
   call require(info == 0, 'upstream lognormal status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_lognormal_bits(i), 'upstream lognormal bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%gumbel(x(1:17), 0.2_dp, 1.3_dp, info)
   call require(info == 0, 'upstream Gumbel status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_gumbel_bits(i), 'upstream Gumbel bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%pareto(x(1:17), 1.7_dp, 2.3_dp, info)
   call require(info == 0, 'upstream Pareto status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_pareto_bits(i), 'upstream Pareto bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%weibull(x(1:17), 1.7_dp, 2.1_dp, info)
   call require(info == 0, 'upstream Weibull status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_weibull_bits(i), 'upstream Weibull bit vector')
   end do
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%skew_normal(x(1:17), 0.4_dp, 1.2_dp, -2.5_dp, info)
   call require(info == 0, 'upstream skew-normal status')
   do i = 1, 17
      call require(transfer(x(i), 0_int64) == expected_skew_normal_bits(i), 'upstream skew-normal bit vector')
   end do

   a = randompack_rng('x256++', seed=777)
   call a%unif(x, -2.0_dp, 3.0_dp, info)
   call require(info == 0 .and. all(x >= -2.0_dp) .and. all(x < 3.0_dp), 'unif')
   call a%normal(x, 2.0_dp, 3.0_dp, info)
   call require(info == 0 .and. all(ieee_is_finite(x)), 'normal')
   call a%skew_normal(x, alpha=2.0_dp, info=info)
   call require(info == 0 .and. all(ieee_is_finite(x)), 'skew normal')
   call a%lognormal(x, info=info)
   call require(info == 0 .and. all(x > 0.0_dp), 'lognormal')
   call a%gumbel(x, info=info)
   call require(info == 0 .and. all(ieee_is_finite(x)), 'gumbel')
   call a%pareto(x, 2.0_dp, 3.0_dp, info)
   call require(info == 0 .and. all(x >= 2.0_dp), 'pareto')
   call a%exp(x, 2.0_dp, info)
   call require(info == 0 .and. all(x >= 0.0_dp), 'exponential')
   call a%gamma(x, 0.7_dp, 2.0_dp, info)
   call require(info == 0 .and. all(x > 0.0_dp), 'gamma')
   call a%chi2(x, 5.0_dp, info)
   call require(info == 0 .and. all(x > 0.0_dp), 'chi2')
   call a%beta(x, 2.0_dp, 5.0_dp, info)
   call require(info == 0 .and. all(x > 0.0_dp) .and. all(x < 1.0_dp), 'beta')
   call a%t(x, 6.0_dp, info)
   call require(info == 0 .and. all(ieee_is_finite(x)), 'student t')
   call a%f(x, 5.0_dp, 8.0_dp, info)
   call require(info == 0 .and. all(x > 0.0_dp), 'F distribution')
   call a%weibull(x, 1.5_dp, 2.0_dp, info)
   call require(info == 0 .and. all(x >= 0.0_dp), 'weibull')

   call a%int(ints, -5, 17, info)
   call require(info == 0 .and. all(ints >= -5) .and. all(ints <= 17), 'integer range')
   call a%perm(p, info)
   call require(info == 0 .and. all([(count(p == i) == 1, i=1,size(p))]), 'permutation')
   call a%sample(40, s, info)
   call require(info == 0 .and. all(s >= 1) .and. all(s <= 40), 'sample bounds')
   do i = 1, size(s)
      call require(count(s == s(i)) == 1, 'sample uniqueness')
   end do
   call a%raw(bytes, info)
   call require(info == 0, 'raw bytes')

   a = randompack_rng('pcg64', seed=123)
   call a%int(ints, -5, 17, info)
   call require(info == 0 .and. all(ints == expected_int), 'upstream bounded integer vector')
   call a%perm(p, info)
   call require(info == 0 .and. all(p == expected_perm), 'upstream permutation vector')
   call a%sample(40, s, info)
   call require(info == 0 .and. all(s == expected_sample), 'upstream sample vector')
   call a%raw(bytes, info)
   call require(info == 0 .and. all(bytes == expected_raw), 'upstream post-discrete raw vector')

   sigma = reshape([1.0_dp, 0.25_dp, 0.25_dp, 2.0_dp], [2, 2])
   mu = [1.0_dp, -1.0_dp]
   call a%mvn(mv, sigma, mu, info)
   call require(info == 0 .and. all(ieee_is_finite(mv)), 'mvn')

   sigma3 = reshape([1.0_dp, 0.25_dp, 0.10_dp, &
      0.25_dp, 2.0_dp, 0.30_dp, 0.10_dp, 0.30_dp, 1.5_dp], [3, 3])
   mu3 = [1.0_dp, -1.0_dp, 0.5_dp]
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%mvn(mv3, sigma3, mu3, info)
   call require(info == 0, 'upstream full-rank MVN status')
   k = 0
   do j = 1, 3
      do i = 1, 5
         k = k + 1
         call require(transfer(mv3(i, j), 0_int64) == expected_mvn_full_bits(k), &
            'upstream full-rank MVN bit vector')
      end do
   end do
   call a%unif(post6, info=info)
   call require(info == 0, 'upstream full-rank post-MVN uniform status')
   do i = 1, 6
      call require(transfer(post6(i), 0_int64) == expected_mvn_full_post_bits(i), &
         'upstream full-rank post-MVN stream')
   end do

   sigma3 = reshape([1.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 4.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, 5.0_dp], [3, 3])
   mu3 = [0.25_dp, -0.5_dp, 1.0_dp]
   a = randompack_rng('pcg64', seed=123, bitexact=.true.)
   call a%mvn(mv3, sigma3, mu3, info)
   call require(info == 0, 'upstream pivoted singular MVN status')
   k = 0
   do j = 1, 3
      do i = 1, 5
         k = k + 1
         call require(transfer(mv3(i, j), 0_int64) == expected_mvn_singular_bits(k), &
            'upstream pivoted singular MVN bit vector')
      end do
   end do
   call a%unif(post6, info=info)
   call require(info == 0, 'upstream singular post-MVN uniform status')
   do i = 1, 6
      call require(transfer(post6(i), 0_int64) == expected_mvn_singular_post_bits(i), &
         'upstream singular post-MVN stream')
   end do

   a = randompack_rng('pcg64', seed=1)
   call a%pcg64_set_inc([3_int64, 0_int64], info)
   call require(info == 0, 'pcg setter')
   a = randompack_rng('cwg128', seed=1)
   call a%cwg128_set_weyl([5_int64, 0_int64], info)
   call require(info == 0, 'cwg setter')
   a = randompack_rng('sfc64', seed=1)
   call a%sfc64_set_abc([1_int64, 2_int64, 3_int64], info)
   call require(info == 0, 'sfc setter')
   a = randompack_rng('chacha20', seed=1)
   call a%chacha_set_nonce([1_int64, 2_int64, 3_int64], info)
   call require(info == 0, 'chacha setter')
   a = randompack_rng('philox', seed=1)
   call a%philox_set_key([1_int64, 2_int64], info)
   call require(info == 0, 'philox setter')
   a = randompack_rng('squares', seed=1)
   call a%squares_set_key(12345_int64, info)
   call require(info == 0, 'squares setter')

   print '(a)', 'All randompack tests passed.'
contains
   subroutine require(condition, message)
      logical, intent(in) :: condition !! Test condition that must evaluate true.
      character(len=*), intent(in) :: message !! Short diagnostic identifying the failed assertion.
      if (.not. condition) then
         write(*, '(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine require
end program test_randompack
