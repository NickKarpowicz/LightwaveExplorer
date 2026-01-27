deviceFunction static inline constexpr float hermite(const float x, const uint8_t n){
    switch(n){
		case 0u:
			return 1;
		case 1u:
			return 2*x;
		case 2u:
			return 4*DeviceFPLib::pow(x, 2) - 2;
		case 3u:
			return 8*DeviceFPLib::pow(x, 3) - 12*x;
		case 4u:
			return 16*DeviceFPLib::pow(x, 4) - 48*DeviceFPLib::pow(x, 2) + 12;
		case 5u:
			return 32*DeviceFPLib::pow(x, 5) - 160*DeviceFPLib::pow(x, 3) + 120*x;
		case 6u:
			return 64*DeviceFPLib::pow(x, 6) - 480*DeviceFPLib::pow(x, 4) + 720*DeviceFPLib::pow(x, 2) - 120;
		case 7u:
			return 128*DeviceFPLib::pow(x, 7) - 1344*DeviceFPLib::pow(x, 5) + 3360*DeviceFPLib::pow(x, 3) - 1680*x;
		case 8u:
			return 256*DeviceFPLib::pow(x, 8) - 3584*DeviceFPLib::pow(x, 6) + 13440*DeviceFPLib::pow(x, 4) - 13440*DeviceFPLib::pow(x, 2) + 1680;
		case 9u:
			return 512*DeviceFPLib::pow(x, 9) - 9216*DeviceFPLib::pow(x, 7) + 48384*DeviceFPLib::pow(x, 5) - 80640*DeviceFPLib::pow(x, 3) + 30240*x;
		case 10u:
			return 1024*DeviceFPLib::pow(x, 10) - 23040*DeviceFPLib::pow(x, 8) + 161280*DeviceFPLib::pow(x, 6) - 403200*DeviceFPLib::pow(x, 4) + 302400*DeviceFPLib::pow(x, 2) - 30240;
		case 11u:
			return 2048*DeviceFPLib::pow(x, 11) - 56320*DeviceFPLib::pow(x, 9) + 506880*DeviceFPLib::pow(x, 7) - 1774080*DeviceFPLib::pow(x, 5) + 2217600*DeviceFPLib::pow(x, 3) - 665280*x;
		case 12u:
			return 4096*DeviceFPLib::pow(x, 12) - 135168*DeviceFPLib::pow(x, 10) + 1520640*DeviceFPLib::pow(x, 8) - 7096320*DeviceFPLib::pow(x, 6) + 13305600*DeviceFPLib::pow(x, 4) - 7983360*DeviceFPLib::pow(x, 2) + 665280;
		case 13u:
			return 8192*DeviceFPLib::pow(x, 13) - 319488*DeviceFPLib::pow(x, 11) + 4392960*DeviceFPLib::pow(x, 9) - 26357760*DeviceFPLib::pow(x, 7) + 69189120*DeviceFPLib::pow(x, 5) - 69189120*DeviceFPLib::pow(x, 3) + 17297280*x;
		case 14u:
			return 16384*DeviceFPLib::pow(x, 14) - 745472*DeviceFPLib::pow(x, 12) + 12300288*DeviceFPLib::pow(x, 10) - 92252160*DeviceFPLib::pow(x, 8) + 322882560*DeviceFPLib::pow(x, 6) - 484323840*DeviceFPLib::pow(x, 4) + 242161920*DeviceFPLib::pow(x, 2) - 17297280;
		case 15u:
			return 32768*DeviceFPLib::pow(x, 15) - 1720320*DeviceFPLib::pow(x, 13) + 33546240*DeviceFPLib::pow(x, 11) - 307507200*DeviceFPLib::pow(x, 9) + 1383782400*DeviceFPLib::pow(x, 7) - 2905943040*DeviceFPLib::pow(x, 5) + 2421619200*DeviceFPLib::pow(x, 3) - 518918400*x;
		case 16u:
			return 65536*DeviceFPLib::pow(x, 16) - 3932160*DeviceFPLib::pow(x, 14) + 89456640*DeviceFPLib::pow(x, 12) - 984023040*DeviceFPLib::pow(x, 10) + 5535129600*DeviceFPLib::pow(x, 8) - 15498362880*DeviceFPLib::pow(x, 6) + 19372953600*DeviceFPLib::pow(x, 4) - 8302694400*DeviceFPLib::pow(x, 2) + 518918400;
		case 17u:
			return 131072*DeviceFPLib::pow(x, 17) - 8912896*DeviceFPLib::pow(x, 15) + 233963520*DeviceFPLib::pow(x, 13) - 3041525760*DeviceFPLib::pow(x, 11) + 20910489600*DeviceFPLib::pow(x, 9) - 75277762560*DeviceFPLib::pow(x, 7) + 131736084480*DeviceFPLib::pow(x, 5) - 94097203200*DeviceFPLib::pow(x, 3) + 17643225600*x;
		case 18u:
			return 262144*DeviceFPLib::pow(x, 18) - 20054016*DeviceFPLib::pow(x, 16) + 601620480*DeviceFPLib::pow(x, 14) - 9124577280*DeviceFPLib::pow(x, 12) + 75277762560*DeviceFPLib::pow(x, 10) - 338749931520*DeviceFPLib::pow(x, 8) + 790416506880*DeviceFPLib::pow(x, 6) - 846874828800*DeviceFPLib::pow(x, 4) + 317578060800*DeviceFPLib::pow(x, 2) - 17643225600;
		case 19u:
			return 524288*DeviceFPLib::pow(x, 19) - 44826624*DeviceFPLib::pow(x, 17) + 1524105216*DeviceFPLib::pow(x, 15) - 26671841280*DeviceFPLib::pow(x, 13) + 260050452480*DeviceFPLib::pow(x, 11) - 1430277488640*DeviceFPLib::pow(x, 9) + 4290832465920*DeviceFPLib::pow(x, 7) - 6436248698880*DeviceFPLib::pow(x, 5) + 4022655436800*DeviceFPLib::pow(x, 3) - 670442572800*x;
		case 20u:
			return 1048576*DeviceFPLib::pow(x, 20) - 99614720*DeviceFPLib::pow(x, 18) + 3810263040*DeviceFPLib::pow(x, 16) - 76205260800*DeviceFPLib::pow(x, 14) + 866834841600*DeviceFPLib::pow(x, 12) - 5721109954560*DeviceFPLib::pow(x, 10) + 21454162329600*DeviceFPLib::pow(x, 8) - 42908324659200*DeviceFPLib::pow(x, 6) + 40226554368000*DeviceFPLib::pow(x, 4) - 13408851456000*DeviceFPLib::pow(x, 2) + 670442572800;
		case 21u:
			return 2097152*DeviceFPLib::pow(x, 21) - 220200960*DeviceFPLib::pow(x, 19) + 9413591040*DeviceFPLib::pow(x, 17) - 213374730240*DeviceFPLib::pow(x, 15) + 2800543334400*DeviceFPLib::pow(x, 13) - 21844238008320*DeviceFPLib::pow(x, 11) + 100119424204800*DeviceFPLib::pow(x, 9) - 257449947955200*DeviceFPLib::pow(x, 7) + 337903056691200*DeviceFPLib::pow(x, 5) - 187723920384000*DeviceFPLib::pow(x, 3) + 28158588057600*x;
		case 22u:
			return 4194304*DeviceFPLib::pow(x, 22) - 484442112*DeviceFPLib::pow(x, 20) + 23011000320*DeviceFPLib::pow(x, 18) - 586780508160*DeviceFPLib::pow(x, 16) + 8801707622400*DeviceFPLib::pow(x, 14) - 80095539363840*DeviceFPLib::pow(x, 12) + 440525466501120*DeviceFPLib::pow(x, 10) - 1415974713753600*DeviceFPLib::pow(x, 8) + 2477955749068800*DeviceFPLib::pow(x, 6) - 2064963124224000*DeviceFPLib::pow(x, 4) + 619488937267200*DeviceFPLib::pow(x, 2) - 28158588057600;
		case 23u:
			return 8388608*DeviceFPLib::pow(x, 23) - 1061158912*DeviceFPLib::pow(x, 21) + 55710842880*DeviceFPLib::pow(x, 19) - 1587759022080*DeviceFPLib::pow(x, 17) + 26991903375360*DeviceFPLib::pow(x, 15) - 283414985441280*DeviceFPLib::pow(x, 13) + 1842197405368320*DeviceFPLib::pow(x, 11) - 7237204092518400*DeviceFPLib::pow(x, 9) + 16283709208166400*DeviceFPLib::pow(x, 7) - 18997660742860800*DeviceFPLib::pow(x, 5) + 9498830371430400*DeviceFPLib::pow(x, 3) - 1295295050649600*x;
		case 24u:
			return 16777216*DeviceFPLib::pow(x, 24) - 2315255808*DeviceFPLib::pow(x, 22) + 133706022912*DeviceFPLib::pow(x, 20) - 4234024058880*DeviceFPLib::pow(x, 18) + 80975710126080*DeviceFPLib::pow(x, 16) - 971708521512960*DeviceFPLib::pow(x, 14) + 7368789621473280*DeviceFPLib::pow(x, 12) - 34738579644088320*DeviceFPLib::pow(x, 10) + 97702255248998400*DeviceFPLib::pow(x, 8) - 151981285942886400*DeviceFPLib::pow(x, 6) + 113985964457164800*DeviceFPLib::pow(x, 4) - 31087081215590400*DeviceFPLib::pow(x, 2) + 1295295050649600;
		case 25u:
			return 33554432*DeviceFPLib::pow(x, 25) - 5033164800*DeviceFPLib::pow(x, 23) + 318347673600*DeviceFPLib::pow(x, 21) - 11142168576000*DeviceFPLib::pow(x, 19) + 238163853312000*DeviceFPLib::pow(x, 17) - 3239028405043200*DeviceFPLib::pow(x, 15) + 28341498544128000*DeviceFPLib::pow(x, 13) - 157902634745856000*DeviceFPLib::pow(x, 11) + 542790306938880000*DeviceFPLib::pow(x, 9) - 1085580613877760000*DeviceFPLib::pow(x, 7) + 1139859644571648000*DeviceFPLib::pow(x, 5) - 518118020259840000*DeviceFPLib::pow(x, 3) + 64764752532480000*x;
		case 26u:
			return 67108864*DeviceFPLib::pow(x, 26) - 10905190400*DeviceFPLib::pow(x, 24) + 752458137600*DeviceFPLib::pow(x, 22) - 28969638297600*DeviceFPLib::pow(x, 20) + 688028909568000*DeviceFPLib::pow(x, 18) - 10526842316390400*DeviceFPLib::pow(x, 16) + 105268423163904000*DeviceFPLib::pow(x, 14) - 684244750565376000*DeviceFPLib::pow(x, 12) + 2822509596082176000*DeviceFPLib::pow(x, 10) - 7056273990205440000*DeviceFPLib::pow(x, 8) + 9878783586287616000*DeviceFPLib::pow(x, 6) - 6735534263377920000*DeviceFPLib::pow(x, 4) + 1683883565844480000*DeviceFPLib::pow(x, 2) - 64764752532480000;
		case 27u:
			return 134217728*DeviceFPLib::pow(x, 27) - 23555211264*DeviceFPLib::pow(x, 25) + 1766640844800*DeviceFPLib::pow(x, 23) - 74493355622400*DeviceFPLib::pow(x, 21) + 1955450585088000*DeviceFPLib::pow(x, 19) - 33438205005004800*DeviceFPLib::pow(x, 17) + 378966323390054400*DeviceFPLib::pow(x, 15) - 2842247425425408000*DeviceFPLib::pow(x, 13) + 13855956198948864000*DeviceFPLib::pow(x, 11) - 42337643941232640000*DeviceFPLib::pow(x, 9) + 76207759094218752000*DeviceFPLib::pow(x, 7) - 72743770044481536000*DeviceFPLib::pow(x, 5) + 30309904185200640000*DeviceFPLib::pow(x, 3) - 3497296636753920000*x;
		case 28u:
			return 268435456*DeviceFPLib::pow(x, 28) - 50734301184*DeviceFPLib::pow(x, 26) + 4122161971200*DeviceFPLib::pow(x, 24) - 189619450675200*DeviceFPLib::pow(x, 22) + 5475261638246400*DeviceFPLib::pow(x, 20) - 104029971126681600*DeviceFPLib::pow(x, 18) + 1326382131865190400*DeviceFPLib::pow(x, 16) - 11368989701701632000*DeviceFPLib::pow(x, 14) + 64661128928428032000*DeviceFPLib::pow(x, 12) - 237090806070902784000*DeviceFPLib::pow(x, 10) + 533454313659531264000*DeviceFPLib::pow(x, 8) - 678941853748494336000*DeviceFPLib::pow(x, 6) + 424338658592808960000*DeviceFPLib::pow(x, 4) - 97924305829109760000*DeviceFPLib::pow(x, 2) + 3497296636753920000;
		case 29u:
			return 536870912*DeviceFPLib::pow(x, 29) - 108984795136*DeviceFPLib::pow(x, 27) + 9563415773184*DeviceFPLib::pow(x, 25) - 478170788659200*DeviceFPLib::pow(x, 23) + 15122151191347200*DeviceFPLib::pow(x, 21) - 317565175018291200*DeviceFPLib::pow(x, 19) + 4525303744010649600*DeviceFPLib::pow(x, 17) - 43960093513246310400*DeviceFPLib::pow(x, 15) + 288488113680678912000*DeviceFPLib::pow(x, 13) - 1250115159282941952000*DeviceFPLib::pow(x, 11) + 3437816688028090368000*DeviceFPLib::pow(x, 9) - 5625518216773238784000*DeviceFPLib::pow(x, 7) + 4922328439676583936000*DeviceFPLib::pow(x, 5) - 1893203246029455360000*DeviceFPLib::pow(x, 3) + 202843204931727360000*x;
		case 30u:
			return 1073741824*DeviceFPLib::pow(x, 30) - 233538846720*DeviceFPLib::pow(x, 28) + 22069421015040*DeviceFPLib::pow(x, 26) - 1195426971648000*DeviceFPLib::pow(x, 24) + 41242230521856000*DeviceFPLib::pow(x, 22) - 952695525054873600*DeviceFPLib::pow(x, 20) + 15084345813368832000*DeviceFPLib::pow(x, 18) - 164850350674673664000*DeviceFPLib::pow(x, 16) + 1236377630060052480000*DeviceFPLib::pow(x, 14) - 6250575796414709760000*DeviceFPLib::pow(x, 12) + 20626900128168542208000*DeviceFPLib::pow(x, 10) - 42191386625799290880000*DeviceFPLib::pow(x, 8) + 49223284396765839360000*DeviceFPLib::pow(x, 6) - 28398048690441830400000*DeviceFPLib::pow(x, 4) + 6085296147951820800000*DeviceFPLib::pow(x, 2) - 202843204931727360000;
		case 31u:
			return 2147483648*DeviceFPLib::pow(x, 31) - 499289948160*DeviceFPLib::pow(x, 29) + 50677929738240*DeviceFPLib::pow(x, 27) - 2964658889687040*DeviceFPLib::pow(x, 25) + 111174708363264000*DeviceFPLib::pow(x, 23) - 2812720121590579200*DeviceFPLib::pow(x, 21) + 49222602127835136000*DeviceFPLib::pow(x, 19) - 601218925989986304000*DeviceFPLib::pow(x, 17) + 5110360870914883584000*DeviceFPLib::pow(x, 15) - 29810438413670154240000*DeviceFPLib::pow(x, 13) + 116260709813313601536000*DeviceFPLib::pow(x, 11) - 290651774533284003840000*DeviceFPLib::pow(x, 9) + 435977661799926005760000*DeviceFPLib::pow(x, 7) - 352135803761478696960000*DeviceFPLib::pow(x, 5) + 125762787057670963200000*DeviceFPLib::pow(x, 3) - 12576278705767096320000*x;
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
			return 4*DeviceFPLib::pow(x, 2) - 2;
		case 3u:
			return 8*DeviceFPLib::pow(x, 3) - 12*x;
		case 4u:
			return 16*DeviceFPLib::pow(x, 4) - 48*DeviceFPLib::pow(x, 2) + 12;
		case 5u:
			return 32*DeviceFPLib::pow(x, 5) - 160*DeviceFPLib::pow(x, 3) + 120*x;
		case 6u:
			return 64*DeviceFPLib::pow(x, 6) - 480*DeviceFPLib::pow(x, 4) + 720*DeviceFPLib::pow(x, 2) - 120;
		case 7u:
			return 128*DeviceFPLib::pow(x, 7) - 1344*DeviceFPLib::pow(x, 5) + 3360*DeviceFPLib::pow(x, 3) - 1680*x;
		case 8u:
			return 256*DeviceFPLib::pow(x, 8) - 3584*DeviceFPLib::pow(x, 6) + 13440*DeviceFPLib::pow(x, 4) - 13440*DeviceFPLib::pow(x, 2) + 1680;
		case 9u:
			return 512*DeviceFPLib::pow(x, 9) - 9216*DeviceFPLib::pow(x, 7) + 48384*DeviceFPLib::pow(x, 5) - 80640*DeviceFPLib::pow(x, 3) + 30240*x;
		case 10u:
			return 1024*DeviceFPLib::pow(x, 10) - 23040*DeviceFPLib::pow(x, 8) + 161280*DeviceFPLib::pow(x, 6) - 403200*DeviceFPLib::pow(x, 4) + 302400*DeviceFPLib::pow(x, 2) - 30240;
		case 11u:
			return 2048*DeviceFPLib::pow(x, 11) - 56320*DeviceFPLib::pow(x, 9) + 506880*DeviceFPLib::pow(x, 7) - 1774080*DeviceFPLib::pow(x, 5) + 2217600*DeviceFPLib::pow(x, 3) - 665280*x;
		case 12u:
			return 4096*DeviceFPLib::pow(x, 12) - 135168*DeviceFPLib::pow(x, 10) + 1520640*DeviceFPLib::pow(x, 8) - 7096320*DeviceFPLib::pow(x, 6) + 13305600*DeviceFPLib::pow(x, 4) - 7983360*DeviceFPLib::pow(x, 2) + 665280;
		case 13u:
			return 8192*DeviceFPLib::pow(x, 13) - 319488*DeviceFPLib::pow(x, 11) + 4392960*DeviceFPLib::pow(x, 9) - 26357760*DeviceFPLib::pow(x, 7) + 69189120*DeviceFPLib::pow(x, 5) - 69189120*DeviceFPLib::pow(x, 3) + 17297280*x;
		case 14u:
			return 16384*DeviceFPLib::pow(x, 14) - 745472*DeviceFPLib::pow(x, 12) + 12300288*DeviceFPLib::pow(x, 10) - 92252160*DeviceFPLib::pow(x, 8) + 322882560*DeviceFPLib::pow(x, 6) - 484323840*DeviceFPLib::pow(x, 4) + 242161920*DeviceFPLib::pow(x, 2) - 17297280;
		case 15u:
			return 32768*DeviceFPLib::pow(x, 15) - 1720320*DeviceFPLib::pow(x, 13) + 33546240*DeviceFPLib::pow(x, 11) - 307507200*DeviceFPLib::pow(x, 9) + 1383782400*DeviceFPLib::pow(x, 7) - 2905943040*DeviceFPLib::pow(x, 5) + 2421619200*DeviceFPLib::pow(x, 3) - 518918400*x;
		case 16u:
			return 65536*DeviceFPLib::pow(x, 16) - 3932160*DeviceFPLib::pow(x, 14) + 89456640*DeviceFPLib::pow(x, 12) - 984023040*DeviceFPLib::pow(x, 10) + 5535129600*DeviceFPLib::pow(x, 8) - 15498362880*DeviceFPLib::pow(x, 6) + 19372953600*DeviceFPLib::pow(x, 4) - 8302694400*DeviceFPLib::pow(x, 2) + 518918400;
		case 17u:
			return 131072*DeviceFPLib::pow(x, 17) - 8912896*DeviceFPLib::pow(x, 15) + 233963520*DeviceFPLib::pow(x, 13) - 3041525760*DeviceFPLib::pow(x, 11) + 20910489600*DeviceFPLib::pow(x, 9) - 75277762560*DeviceFPLib::pow(x, 7) + 131736084480*DeviceFPLib::pow(x, 5) - 94097203200*DeviceFPLib::pow(x, 3) + 17643225600*x;
		case 18u:
			return 262144*DeviceFPLib::pow(x, 18) - 20054016*DeviceFPLib::pow(x, 16) + 601620480*DeviceFPLib::pow(x, 14) - 9124577280*DeviceFPLib::pow(x, 12) + 75277762560*DeviceFPLib::pow(x, 10) - 338749931520*DeviceFPLib::pow(x, 8) + 790416506880*DeviceFPLib::pow(x, 6) - 846874828800*DeviceFPLib::pow(x, 4) + 317578060800*DeviceFPLib::pow(x, 2) - 17643225600;
		case 19u:
			return 524288*DeviceFPLib::pow(x, 19) - 44826624*DeviceFPLib::pow(x, 17) + 1524105216*DeviceFPLib::pow(x, 15) - 26671841280*DeviceFPLib::pow(x, 13) + 260050452480*DeviceFPLib::pow(x, 11) - 1430277488640*DeviceFPLib::pow(x, 9) + 4290832465920*DeviceFPLib::pow(x, 7) - 6436248698880*DeviceFPLib::pow(x, 5) + 4022655436800*DeviceFPLib::pow(x, 3) - 670442572800*x;
		case 20u:
			return 1048576*DeviceFPLib::pow(x, 20) - 99614720*DeviceFPLib::pow(x, 18) + 3810263040*DeviceFPLib::pow(x, 16) - 76205260800*DeviceFPLib::pow(x, 14) + 866834841600*DeviceFPLib::pow(x, 12) - 5721109954560*DeviceFPLib::pow(x, 10) + 21454162329600*DeviceFPLib::pow(x, 8) - 42908324659200*DeviceFPLib::pow(x, 6) + 40226554368000*DeviceFPLib::pow(x, 4) - 13408851456000*DeviceFPLib::pow(x, 2) + 670442572800;
		case 21u:
			return 2097152*DeviceFPLib::pow(x, 21) - 220200960*DeviceFPLib::pow(x, 19) + 9413591040*DeviceFPLib::pow(x, 17) - 213374730240*DeviceFPLib::pow(x, 15) + 2800543334400*DeviceFPLib::pow(x, 13) - 21844238008320*DeviceFPLib::pow(x, 11) + 100119424204800*DeviceFPLib::pow(x, 9) - 257449947955200*DeviceFPLib::pow(x, 7) + 337903056691200*DeviceFPLib::pow(x, 5) - 187723920384000*DeviceFPLib::pow(x, 3) + 28158588057600*x;
		case 22u:
			return 4194304*DeviceFPLib::pow(x, 22) - 484442112*DeviceFPLib::pow(x, 20) + 23011000320*DeviceFPLib::pow(x, 18) - 586780508160*DeviceFPLib::pow(x, 16) + 8801707622400*DeviceFPLib::pow(x, 14) - 80095539363840*DeviceFPLib::pow(x, 12) + 440525466501120*DeviceFPLib::pow(x, 10) - 1415974713753600*DeviceFPLib::pow(x, 8) + 2477955749068800*DeviceFPLib::pow(x, 6) - 2064963124224000*DeviceFPLib::pow(x, 4) + 619488937267200*DeviceFPLib::pow(x, 2) - 28158588057600;
		case 23u:
			return 8388608*DeviceFPLib::pow(x, 23) - 1061158912*DeviceFPLib::pow(x, 21) + 55710842880*DeviceFPLib::pow(x, 19) - 1587759022080*DeviceFPLib::pow(x, 17) + 26991903375360*DeviceFPLib::pow(x, 15) - 283414985441280*DeviceFPLib::pow(x, 13) + 1842197405368320*DeviceFPLib::pow(x, 11) - 7237204092518400*DeviceFPLib::pow(x, 9) + 16283709208166400*DeviceFPLib::pow(x, 7) - 18997660742860800*DeviceFPLib::pow(x, 5) + 9498830371430400*DeviceFPLib::pow(x, 3) - 1295295050649600*x;
		case 24u:
			return 16777216*DeviceFPLib::pow(x, 24) - 2315255808*DeviceFPLib::pow(x, 22) + 133706022912*DeviceFPLib::pow(x, 20) - 4234024058880*DeviceFPLib::pow(x, 18) + 80975710126080*DeviceFPLib::pow(x, 16) - 971708521512960*DeviceFPLib::pow(x, 14) + 7368789621473280*DeviceFPLib::pow(x, 12) - 34738579644088320*DeviceFPLib::pow(x, 10) + 97702255248998400*DeviceFPLib::pow(x, 8) - 151981285942886400*DeviceFPLib::pow(x, 6) + 113985964457164800*DeviceFPLib::pow(x, 4) - 31087081215590400*DeviceFPLib::pow(x, 2) + 1295295050649600;
		case 25u:
			return 33554432*DeviceFPLib::pow(x, 25) - 5033164800*DeviceFPLib::pow(x, 23) + 318347673600*DeviceFPLib::pow(x, 21) - 11142168576000*DeviceFPLib::pow(x, 19) + 238163853312000*DeviceFPLib::pow(x, 17) - 3239028405043200*DeviceFPLib::pow(x, 15) + 28341498544128000*DeviceFPLib::pow(x, 13) - 157902634745856000*DeviceFPLib::pow(x, 11) + 542790306938880000*DeviceFPLib::pow(x, 9) - 1085580613877760000*DeviceFPLib::pow(x, 7) + 1139859644571648000*DeviceFPLib::pow(x, 5) - 518118020259840000*DeviceFPLib::pow(x, 3) + 64764752532480000*x;
		case 26u:
			return 67108864*DeviceFPLib::pow(x, 26) - 10905190400*DeviceFPLib::pow(x, 24) + 752458137600*DeviceFPLib::pow(x, 22) - 28969638297600*DeviceFPLib::pow(x, 20) + 688028909568000*DeviceFPLib::pow(x, 18) - 10526842316390400*DeviceFPLib::pow(x, 16) + 105268423163904000*DeviceFPLib::pow(x, 14) - 684244750565376000*DeviceFPLib::pow(x, 12) + 2822509596082176000*DeviceFPLib::pow(x, 10) - 7056273990205440000*DeviceFPLib::pow(x, 8) + 9878783586287616000*DeviceFPLib::pow(x, 6) - 6735534263377920000*DeviceFPLib::pow(x, 4) + 1683883565844480000*DeviceFPLib::pow(x, 2) - 64764752532480000;
		case 27u:
			return 134217728*DeviceFPLib::pow(x, 27) - 23555211264*DeviceFPLib::pow(x, 25) + 1766640844800*DeviceFPLib::pow(x, 23) - 74493355622400*DeviceFPLib::pow(x, 21) + 1955450585088000*DeviceFPLib::pow(x, 19) - 33438205005004800*DeviceFPLib::pow(x, 17) + 378966323390054400*DeviceFPLib::pow(x, 15) - 2842247425425408000*DeviceFPLib::pow(x, 13) + 13855956198948864000*DeviceFPLib::pow(x, 11) - 42337643941232640000*DeviceFPLib::pow(x, 9) + 76207759094218752000*DeviceFPLib::pow(x, 7) - 72743770044481536000*DeviceFPLib::pow(x, 5) + 30309904185200640000*DeviceFPLib::pow(x, 3) - 3497296636753920000*x;
		case 28u:
			return 268435456*DeviceFPLib::pow(x, 28) - 50734301184*DeviceFPLib::pow(x, 26) + 4122161971200*DeviceFPLib::pow(x, 24) - 189619450675200*DeviceFPLib::pow(x, 22) + 5475261638246400*DeviceFPLib::pow(x, 20) - 104029971126681600*DeviceFPLib::pow(x, 18) + 1326382131865190400*DeviceFPLib::pow(x, 16) - 11368989701701632000*DeviceFPLib::pow(x, 14) + 64661128928428032000*DeviceFPLib::pow(x, 12) - 237090806070902784000*DeviceFPLib::pow(x, 10) + 533454313659531264000*DeviceFPLib::pow(x, 8) - 678941853748494336000*DeviceFPLib::pow(x, 6) + 424338658592808960000*DeviceFPLib::pow(x, 4) - 97924305829109760000*DeviceFPLib::pow(x, 2) + 3497296636753920000;
		case 29u:
			return 536870912*DeviceFPLib::pow(x, 29) - 108984795136*DeviceFPLib::pow(x, 27) + 9563415773184*DeviceFPLib::pow(x, 25) - 478170788659200*DeviceFPLib::pow(x, 23) + 15122151191347200*DeviceFPLib::pow(x, 21) - 317565175018291200*DeviceFPLib::pow(x, 19) + 4525303744010649600*DeviceFPLib::pow(x, 17) - 43960093513246310400*DeviceFPLib::pow(x, 15) + 288488113680678912000*DeviceFPLib::pow(x, 13) - 1250115159282941952000*DeviceFPLib::pow(x, 11) + 3437816688028090368000*DeviceFPLib::pow(x, 9) - 5625518216773238784000*DeviceFPLib::pow(x, 7) + 4922328439676583936000*DeviceFPLib::pow(x, 5) - 1893203246029455360000*DeviceFPLib::pow(x, 3) + 202843204931727360000*x;
		case 30u:
			return 1073741824*DeviceFPLib::pow(x, 30) - 233538846720*DeviceFPLib::pow(x, 28) + 22069421015040*DeviceFPLib::pow(x, 26) - 1195426971648000*DeviceFPLib::pow(x, 24) + 41242230521856000*DeviceFPLib::pow(x, 22) - 952695525054873600*DeviceFPLib::pow(x, 20) + 15084345813368832000*DeviceFPLib::pow(x, 18) - 164850350674673664000*DeviceFPLib::pow(x, 16) + 1236377630060052480000*DeviceFPLib::pow(x, 14) - 6250575796414709760000*DeviceFPLib::pow(x, 12) + 20626900128168542208000*DeviceFPLib::pow(x, 10) - 42191386625799290880000*DeviceFPLib::pow(x, 8) + 49223284396765839360000*DeviceFPLib::pow(x, 6) - 28398048690441830400000*DeviceFPLib::pow(x, 4) + 6085296147951820800000*DeviceFPLib::pow(x, 2) - 202843204931727360000;
		case 31u:
			return 2147483648*DeviceFPLib::pow(x, 31) - 499289948160*DeviceFPLib::pow(x, 29) + 50677929738240*DeviceFPLib::pow(x, 27) - 2964658889687040*DeviceFPLib::pow(x, 25) + 111174708363264000*DeviceFPLib::pow(x, 23) - 2812720121590579200*DeviceFPLib::pow(x, 21) + 49222602127835136000*DeviceFPLib::pow(x, 19) - 601218925989986304000*DeviceFPLib::pow(x, 17) + 5110360870914883584000*DeviceFPLib::pow(x, 15) - 29810438413670154240000*DeviceFPLib::pow(x, 13) + 116260709813313601536000*DeviceFPLib::pow(x, 11) - 290651774533284003840000*DeviceFPLib::pow(x, 9) + 435977661799926005760000*DeviceFPLib::pow(x, 7) - 352135803761478696960000*DeviceFPLib::pow(x, 5) + 125762787057670963200000*DeviceFPLib::pow(x, 3) - 12576278705767096320000*x;
		default: return 0.0;
    }
}
