deviceFunction static inline constexpr float hermite(const float x, const uint8_t n){
    switch(n){
		case 0u:
			return 1;
		case 1u:
			return 2*x;
		case 2u:
			return 4*deviceFPLib::pow(x, 2) - 2;
		case 3u:
			return 8*deviceFPLib::pow(x, 3) - 12*x;
		case 4u:
			return 16*deviceFPLib::pow(x, 4) - 48*deviceFPLib::pow(x, 2) + 12;
		case 5u:
			return 32*deviceFPLib::pow(x, 5) - 160*deviceFPLib::pow(x, 3) + 120*x;
		case 6u:
			return 64*deviceFPLib::pow(x, 6) - 480*deviceFPLib::pow(x, 4) + 720*deviceFPLib::pow(x, 2) - 120;
		case 7u:
			return 128*deviceFPLib::pow(x, 7) - 1344*deviceFPLib::pow(x, 5) + 3360*deviceFPLib::pow(x, 3) - 1680*x;
		case 8u:
			return 256*deviceFPLib::pow(x, 8) - 3584*deviceFPLib::pow(x, 6) + 13440*deviceFPLib::pow(x, 4) - 13440*deviceFPLib::pow(x, 2) + 1680;
		case 9u:
			return 512*deviceFPLib::pow(x, 9) - 9216*deviceFPLib::pow(x, 7) + 48384*deviceFPLib::pow(x, 5) - 80640*deviceFPLib::pow(x, 3) + 30240*x;
		case 10u:
			return 1024*deviceFPLib::pow(x, 10) - 23040*deviceFPLib::pow(x, 8) + 161280*deviceFPLib::pow(x, 6) - 403200*deviceFPLib::pow(x, 4) + 302400*deviceFPLib::pow(x, 2) - 30240;
		case 11u:
			return 2048*deviceFPLib::pow(x, 11) - 56320*deviceFPLib::pow(x, 9) + 506880*deviceFPLib::pow(x, 7) - 1774080*deviceFPLib::pow(x, 5) + 2217600*deviceFPLib::pow(x, 3) - 665280*x;
		case 12u:
			return 4096*deviceFPLib::pow(x, 12) - 135168*deviceFPLib::pow(x, 10) + 1520640*deviceFPLib::pow(x, 8) - 7096320*deviceFPLib::pow(x, 6) + 13305600*deviceFPLib::pow(x, 4) - 7983360*deviceFPLib::pow(x, 2) + 665280;
		case 13u:
			return 8192*deviceFPLib::pow(x, 13) - 319488*deviceFPLib::pow(x, 11) + 4392960*deviceFPLib::pow(x, 9) - 26357760*deviceFPLib::pow(x, 7) + 69189120*deviceFPLib::pow(x, 5) - 69189120*deviceFPLib::pow(x, 3) + 17297280*x;
		case 14u:
			return 16384*deviceFPLib::pow(x, 14) - 745472*deviceFPLib::pow(x, 12) + 12300288*deviceFPLib::pow(x, 10) - 92252160*deviceFPLib::pow(x, 8) + 322882560*deviceFPLib::pow(x, 6) - 484323840*deviceFPLib::pow(x, 4) + 242161920*deviceFPLib::pow(x, 2) - 17297280;
		case 15u:
			return 32768*deviceFPLib::pow(x, 15) - 1720320*deviceFPLib::pow(x, 13) + 33546240*deviceFPLib::pow(x, 11) - 307507200*deviceFPLib::pow(x, 9) + 1383782400*deviceFPLib::pow(x, 7) - 2905943040*deviceFPLib::pow(x, 5) + 2421619200*deviceFPLib::pow(x, 3) - 518918400*x;
		default: return 0.0f;
    }
}
deviceFunction static inline constexpr double hermite(const double x, const uint8_t n){
    switch(n){
		case 0u:
			return 1;
		case 1u:
			return 2*x;
		case 2u:
			return 4*deviceFPLib::pow(x, 2) - 2;
		case 3u:
			return 8*deviceFPLib::pow(x, 3) - 12*x;
		case 4u:
			return 16*deviceFPLib::pow(x, 4) - 48*deviceFPLib::pow(x, 2) + 12;
		case 5u:
			return 32*deviceFPLib::pow(x, 5) - 160*deviceFPLib::pow(x, 3) + 120*x;
		case 6u:
			return 64*deviceFPLib::pow(x, 6) - 480*deviceFPLib::pow(x, 4) + 720*deviceFPLib::pow(x, 2) - 120;
		case 7u:
			return 128*deviceFPLib::pow(x, 7) - 1344*deviceFPLib::pow(x, 5) + 3360*deviceFPLib::pow(x, 3) - 1680*x;
		case 8u:
			return 256*deviceFPLib::pow(x, 8) - 3584*deviceFPLib::pow(x, 6) + 13440*deviceFPLib::pow(x, 4) - 13440*deviceFPLib::pow(x, 2) + 1680;
		case 9u:
			return 512*deviceFPLib::pow(x, 9) - 9216*deviceFPLib::pow(x, 7) + 48384*deviceFPLib::pow(x, 5) - 80640*deviceFPLib::pow(x, 3) + 30240*x;
		case 10u:
			return 1024*deviceFPLib::pow(x, 10) - 23040*deviceFPLib::pow(x, 8) + 161280*deviceFPLib::pow(x, 6) - 403200*deviceFPLib::pow(x, 4) + 302400*deviceFPLib::pow(x, 2) - 30240;
		case 11u:
			return 2048*deviceFPLib::pow(x, 11) - 56320*deviceFPLib::pow(x, 9) + 506880*deviceFPLib::pow(x, 7) - 1774080*deviceFPLib::pow(x, 5) + 2217600*deviceFPLib::pow(x, 3) - 665280*x;
		case 12u:
			return 4096*deviceFPLib::pow(x, 12) - 135168*deviceFPLib::pow(x, 10) + 1520640*deviceFPLib::pow(x, 8) - 7096320*deviceFPLib::pow(x, 6) + 13305600*deviceFPLib::pow(x, 4) - 7983360*deviceFPLib::pow(x, 2) + 665280;
		case 13u:
			return 8192*deviceFPLib::pow(x, 13) - 319488*deviceFPLib::pow(x, 11) + 4392960*deviceFPLib::pow(x, 9) - 26357760*deviceFPLib::pow(x, 7) + 69189120*deviceFPLib::pow(x, 5) - 69189120*deviceFPLib::pow(x, 3) + 17297280*x;
		case 14u:
			return 16384*deviceFPLib::pow(x, 14) - 745472*deviceFPLib::pow(x, 12) + 12300288*deviceFPLib::pow(x, 10) - 92252160*deviceFPLib::pow(x, 8) + 322882560*deviceFPLib::pow(x, 6) - 484323840*deviceFPLib::pow(x, 4) + 242161920*deviceFPLib::pow(x, 2) - 17297280;
		case 15u:
			return 32768*deviceFPLib::pow(x, 15) - 1720320*deviceFPLib::pow(x, 13) + 33546240*deviceFPLib::pow(x, 11) - 307507200*deviceFPLib::pow(x, 9) + 1383782400*deviceFPLib::pow(x, 7) - 2905943040*deviceFPLib::pow(x, 5) + 2421619200*deviceFPLib::pow(x, 3) - 518918400*x;
		default: return 0.0;
    }
}
