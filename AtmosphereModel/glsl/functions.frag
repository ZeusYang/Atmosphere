/*
The functions provided in this file are organized as follows:
Transmittance
    Computation
    Precomputation
    Lookup
Single scattering
    Computation
    Precomputation
    Lookup
Multiple scattering
    Computation
    Precomputation
    Lookup
Ground irradiance
    Computation
    Precomputation
    Lookup
Rendering
    Sky
    Aerial perspective
    Ground
*/

/**
 * 以下几个clamp函数用以限制数值在其对应的数值域内
 */

Number ClampCosine(Number mu){
	return(clamp(mu, Number(-1.0), Number(1.0)));
}


Length ClampDistance(Length d){
	return(max(d, 0.0 * m));
}


Length ClampRadius( IN(AtmosphereParameters)atmosphere, Length r ){
	return(clamp( r, atmosphere.bottom_radius, atmosphere.top_radius ) );
}


Length SafeSqrt( Area a ){
	return(sqrt( max( a, 0.0 * m2 ) ) );
}

/**
 * 功能:
 *  用求根公式求解二元一次方程x^2+2urx+r^2-t^2=0
 *  其中u(即下面的mu)是视线天顶角的cos值,而t(即下面的top_radius)是大气层最外层半径
 *  r是视点位置向量的z分量,求解的根是视点p沿视线到大气层顶层的距离
 * 传入参数：
 *  atmosphere大气模型参数(其中top_radius要用到),r为视点高度,mu是视线天顶角的cos值
 **/
Length DistanceToTopAtmosphereBoundary( IN(AtmosphereParameters)atmosphere,
					Length r, Number mu ){
	Area discriminant = r*r*(mu*mu-1.0)+atmosphere.top_radius*atmosphere.top_radius;//判别式
	return(ClampDistance( -r * mu + SafeSqrt( discriminant ) ) ); /* 这里是方程解中"+"的那个根 */
}

/**
 * 功能:
 *  与函数DistanceToTopAtmosphereBoundary类似
 *  不过这里计算的是视点沿视线方向到地球表面交点的距离
 * 传入参数：
 *  atmosphere大气模型参数(其中top_radius要用到),r为视点高度,mu是视线天顶角的cos值
 **/
Length DistanceToBottomAtmosphereBoundary( IN(AtmosphereParameters)atmosphere,
					   Length r, Number mu ){
	Area discriminant = r*r*(mu*mu-1.0)+atmosphere.bottom_radius*atmosphere.bottom_radius;
	return(ClampDistance( -r * mu - SafeSqrt( discriminant ) ) ); /* 这里是方程解中"-"的那个根 */
}

/**
 * 功能:
 *  设i为视点沿视线方向与地球表面的交点,视点所在为p
 *  则当二元一次方程d^2+2rud+r^2-t^2=0有解d>=0时,射线pi与地面有交点
 *  该方程意义前面已经说明,求得是距离,而距离一定>=0,故只需该方程有解，那么就有交点
 *  这个函数判断通过方程的判别式判断方程是否有解
 * 传入参数：
 *  atmosphere大气模型参数(其中top_radius要用到),r为视点高度,mu是视线天顶角的cos值
 **/
bool RayIntersectsGround( IN(AtmosphereParameters)atmosphere,
			  Length r, Number mu ){
	// 天顶角大于90度(cos值<0) 且 判别式>=0
	return(mu < 0.0 && r * r * (mu * mu - 1.0) +
	       atmosphere.bottom_radius * atmosphere.bottom_radius >= 0.0 * m2);
}

/**
 * 功能:
 *  根据给定高度计算密度,密度='exp_term'*exp('exp_scale' * h)+'linear_term'*h+'constant_term'
 * 传入参数：
 *  layer大气层,altitude为海拔高度
 **/
Number GetLayerDensity( IN(DensityProfileLayer)layer, Length altitude ){
	/* 密度计算'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term' */
	Number density = layer.exp_term * exp( layer.exp_scale * altitude ) +
			 layer.linear_term * altitude + layer.constant_term;
	return(clamp( density, Number( 0.0 ), Number( 1.0 ) ) );
}

/**
 * 功能:
 *  根据给定海拔高度，计算相应的密度值，整个大气层的密度模型有两层
 * 传入参数：
 *  profile为大气层,altitude为海拔高度
 **/
Number GetProfileDensity( IN(DensityProfile)profile, Length altitude ){
	return(altitude < profile.layers[0].width ?
	       GetLayerDensity( profile.layers[0], altitude ) :
	       GetLayerDensity( profile.layers[1], altitude ) );
}

/**
 * 功能:
 *  用梯度法和光线步进法(Ray Marching)计算从p到大气层顶层点的光学长度
 * 传入参数：
 *  atmosphere为大气模型,profile为相应的密度分布,r为视点高度,mu是视线天顶角的cos值
 **/
Length ComputeOpticalLengthToTopAtmosphereBoundary(
	IN(AtmosphereParameters)atmosphere, IN(DensityProfile)profile,
	Length r, Number mu ){
	/* 取500个采样进行积分 */
	const int SAMPLE_COUNT = 500;
	/* 积分区间步长 */
	Length dx = DistanceToTopAtmosphereBoundary(atmosphere,r,mu)/Number(SAMPLE_COUNT);
	/* 用梯度法循环计算进行积分 */
	Length result = 0.0 * m;
	for ( int i = 0; i <= SAMPLE_COUNT; ++i) {
		Length d_i = Number( i ) * dx;
		/* 当前采样点到星球中心的距离 */
		Length r_i = sqrt( d_i * d_i + 2.0 * r * mu * d_i + r * r );
		/* 获取当前采样点的密度值,注意传入的是海拔高度 */
		Number y_i = GetProfileDensity( profile, r_i - atmosphere.bottom_radius );
		/* (根据梯度计算法则,在积分段两端取0.5的权值，其余为1.0 */
		Number weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
		result += y_i * weight_i * dx;
	}
	return result;
}

/**
 * 功能:
 *  计算视点沿视线方向与大气层顶部交点的透射率(或称光学深度)
 *  结果为三部分叠加取负的exp值:rayleigh光学长度+mie光学长度+吸收光线分子的光学长度
 *  对应的介质分别为:空气分子、气溶胶、吸收光线的介质
 * 传入参数：
 *  atmosphere为大气模型,r为视点高度,mu是视线天顶角的cos值
 **/
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundary(
	IN(AtmosphereParameters)atmosphere, Length r, Number mu ){
	return(exp( -(
			    atmosphere.rayleigh_scattering *
			    ComputeOpticalLengthToTopAtmosphereBoundary(
				    atmosphere, atmosphere.rayleigh_density, r, mu ) +
			    atmosphere.mie_extinction *
			    ComputeOpticalLengthToTopAtmosphereBoundary(
				    atmosphere, atmosphere.mie_density, r, mu ) +
			    atmosphere.absorption_extinction *
			    ComputeOpticalLengthToTopAtmosphereBoundary(
				    atmosphere, atmosphere.absorption_density, r, mu ) ) ) );
}

/**
 * 功能:
 *  将[0,1]的x映射到[0.5/n,1.0-0.5/n],其中n是纹理大小
 *  原因是防止在纹理边界部分采样产生一些外推值
 * 传入参数：
 *  x要映射的值,texture_size纹理大小
 **/
Number GetTextureCoordFromUnitRange( Number x, int texture_size ){
	return(0.5 / Number( texture_size ) + x * (1.0 - 1.0 / Number( texture_size ) ) );
}

/**
 * 功能:
 *  GetTextureCoordFromUnitRange的逆过程
 * 传入参数：
 *  u要映射的值,texture_size纹理大小
 **/
Number GetUnitRangeFromTextureCoord( Number u, int texture_size ){
	return( (u - 0.5 / Number( texture_size ) ) / (1.0 - 1.0 / Number( texture_size ) ) );
}

/**
 * 功能:
 *  将(r,mu)通过哈希函数映射到纹理坐标(u,v)
 * 传入参数：
 *  atmosphere大气模型参数,r为视点高度,mu是视线天顶角的cos值
 **/
vec2 GetTransmittanceTextureUvFromRMu(IN(AtmosphereParameters)atmosphere,
				       Length r, Number mu){
	/* 地表切线射线与大气层顶部的交点与切线的起点之间的距离 */
	Length H = sqrt( atmosphere.top_radius * atmosphere.top_radius -
			 atmosphere.bottom_radius * atmosphere.bottom_radius );
	/* 过视点的地表切线的切点到视点的距离 */
	Length rho = SafeSqrt(r*r-atmosphere.bottom_radius*atmosphere.bottom_radius);
	/* d为视点沿视线方向与大气层外层的交点的距离,后两个是d的上界、下界*/
	Length	d		= DistanceToTopAtmosphereBoundary( atmosphere, r, mu );
	Length	d_min	= atmosphere.top_radius - r;// 下界是视点到大气层的垂直距离
	Length	d_max	= rho + H; // 上界的情况是视线刚好与地表相切
	Number	x_mu	= (d - d_min) / (d_max - d_min); // 将mu映射到[0,1]
	Number	x_r		= rho / H; 
	return vec2(GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH),
		     	  GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT));
}

/**
 * 功能:
 *  将(u,v)映射到纹理坐标(r,mu),即上面函数的逆过程
 * 传入参数：
 *  atmosphere大气模型参数,r为视点高度,mu是视线天顶角的cos值,uv是要变换的纹理坐标
 **/
void GetRMuFromTransmittanceTextureUv( IN(AtmosphereParameters)atmosphere,
				       IN(vec2)uv, OUT(Length)r, OUT(Number)mu ){
	/* 变换回来 */
	Number	x_mu	= GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
	Number	x_r		= GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);
	/* 地表切线射线与大气层顶部的交点与切线的起点之间的距离 */
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
			 		   atmosphere.bottom_radius * atmosphere.bottom_radius);
	/* 过视点的地表切线的切点到视点的距离,用于计算r(三角形勾股定理) */
	Length rho = H * x_r;
	r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
	/* d为视点沿视线方向与大气层外层的交点的距离,这里计算d的上、下界,从而得出d */
	Length	d_min	= atmosphere.top_radius - r;
	Length	d_max	= rho + H;
	Length	d		= d_min + x_mu * (d_max - d_min);
   // 根据H、rho、d推出视线天顶角cos值，由前面提到的求d方程d^2+2rud+r^2-t^2=0推出其中的u
	mu	= d == 0.0 * m ? Number(1.0) : (H * H - rho * rho - d * d) / (2.0 * r * d);
	mu	= ClampCosine(mu);
}

/**
 * 功能:
 *  将传入的纹理坐标(u,v)映射回(r,mu),然后(r,mu)用于计算光学深度
 * 传入参数：
 *  atmosphere大气参数模型,frag_coord为当前的片元坐标
 **/
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundaryTexture(
	IN(AtmosphereParameters)atmosphere, IN(vec2)frag_coord ){
	const vec2 TRANSMITTANCE_TEXTURE_SIZE =
				vec2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
	Length	r;/* r为视点高度,mu是视线天顶角的cos值 */
	Number	mu;
	GetRMuFromTransmittanceTextureUv(atmosphere, 
						frag_coord/TRANSMITTANCE_TEXTURE_SIZE, r, mu );
	return	ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, mu);
}

/**
 * 功能:
 *  将(r,mu)映射到纹理坐标(u,v),然后根据(u,v)获取相应纹理单元存储的光学深度
 *  即用(r,mu)去查找预计算存放的纹理,返回相应单元的值
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture透射率纹理,
 *  r为视点高度,mu是视线方向天顶角的cos值
 **/
DimensionlessSpectrum GetTransmittanceToTopAtmosphereBoundary(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu){
	vec2 uv = GetTransmittanceTextureUvFromRMu(atmosphere, r, mu);
	return DimensionlessSpectrum(texture(transmittance_texture, uv));
}

/**
 * 功能:
 *  计算视点p与从视点触发沿视线长度为d的点q之间的光学深度
 *  pq光学深度=pi光学深度/qi光学深度,i是p沿视线与大气顶层的交点
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture光学深度纹理,r为视点p高度
 *  ,mu是视线天顶角的cos值,d是向量pq的长度,ray_r_mu_intersects_ground射线是否与地面相交
 **/
DimensionlessSpectrum GetTransmittance(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu, Length d, bool ray_r_mu_intersects_ground){
	/* 计算点q处的高度与天顶角cos值 */
	Length	r_d		= ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
	Number	mu_d	= ClampCosine((r * mu + d) / r_d);// oq向量与pi向量点积除以两向量长度之积

	/* 如果射线pi与地面有交点,取反向的天顶角 */
	if ( ray_r_mu_intersects_ground ){
		return min(
			       GetTransmittanceToTopAtmosphereBoundary(
				       atmosphere, transmittance_texture, r_d, -mu_d ) /
			       GetTransmittanceToTopAtmosphereBoundary(
				       atmosphere, transmittance_texture, r, -mu ),
			       DimensionlessSpectrum(1.0));
	} else {
		return min(
			       GetTransmittanceToTopAtmosphereBoundary(
				       atmosphere, transmittance_texture, r, mu ) /
			       GetTransmittanceToTopAtmosphereBoundary(
				       atmosphere, transmittance_texture, r_d, mu_d ),
			       DimensionlessSpectrum(1.0));
	}
}

/**
 * 功能:
 *  计算某一点到太阳的光学
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture透射率纹理,r为视点p海拔高度
 *  ,mu_s是射线与太阳夹角的cos值
 **/
DimensionlessSpectrum GetTransmittanceToSun(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu_s ){
	/* 计算过视点p的地面切线射线的天顶角的sin和cos值 */
	Number	sin_theta_h	= atmosphere.bottom_radius / r;
	Number	cos_theta_h	= -sqrt(max(1.0-sin_theta_h*sin_theta_h, 0.0));//取负是因为theta一定大于等于90度
	/* 到太阳的透射率=光学长度*太阳在地平面以上部分的区域部分 */
	return(GetTransmittanceToTopAtmosphereBoundary(
		       atmosphere, transmittance_texture, r, mu_s ) *
	       smoothstep(-sin_theta_h * atmosphere.sun_angular_radius/rad,// rad为弧度单位
			   			  sin_theta_h * atmosphere.sun_angular_radius/rad,
			   			  mu_s-cos_theta_h));
}

/**
 * 功能:
 *  计算单次散射积分,返回rayleigh散射系数和mie散射系数积分值,即光学深度(不是光学长度)
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture光学长度纹理,r为视点p高度
 *  mu为视线天顶角cos值,mu_s是太阳方向天顶角的cos值,d为视点p沿视线方向与q的距离,
 *  nu是视线与太阳单位方向向量的夹角cos值,ray_r_mu_intersects_ground射线是否与地面相交
 **/
void ComputeSingleScatteringIntegrand(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu, Number mu_s, Number nu, Length d,
	bool ray_r_mu_intersects_ground,
	OUT(DimensionlessSpectrum)rayleigh, OUT(DimensionlessSpectrum)mie ){
	/* 设q为射线pi上的一点,r_d为q的海拔高度 */
	Length r_d 	 = ClampRadius( atmosphere, sqrt( d * d + 2.0 * r * mu * d + r * r ) );
	/* mu_s_d为向量0q与太阳单位方向向量的夹角的cos值 */
	Number mu_s_d = ClampCosine( (r * mu_s + d * nu) / r_d );
	/* 相乘是因为GetTransmittance返回的是exp,相乘即他们的指数相加 */
	DimensionlessSpectrum transmittance =
		GetTransmittance(atmosphere, transmittance_texture, r, mu, d,ray_r_mu_intersects_ground) *
		GetTransmittanceToSun(atmosphere, transmittance_texture, r_d, mu_s_d);
	rayleigh = transmittance * GetProfileDensity(
		atmosphere.rayleigh_density, r_d - atmosphere.bottom_radius);
	mie 	  = transmittance * GetProfileDensity(
		atmosphere.mie_density, r_d - atmosphere.bottom_radius);
}

/**
 * 功能:
 *  计算射线(r,mu)到大气层边界最近的距离
 * 传入参数：
 *  atmosphere大气参数模型,r为视点p高度,mu为视线天顶角cos值,
 *	 ray_r_mu_intersects_ground射线是否与地面相交
 **/
Length DistanceToNearestAtmosphereBoundary( IN(AtmosphereParameters)atmosphere,
					    Length r, Number mu, bool ray_r_mu_intersects_ground ){
	/* 若射线与地面有交点,返回射点到地球表面的距离 */
	if ( ray_r_mu_intersects_ground ){
		return DistanceToBottomAtmosphereBoundary(atmosphere, r, mu);
	}else{ /* 否则返回到大气层顶部的距离 */
		return DistanceToTopAtmosphereBoundary( atmosphere, r, mu);
	}
}

/**
 * 功能:
 *  用梯度法和光线步进计算单次内散射积分
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture光学长度纹理,r为视点p高度
 *  mu为视线天顶角cos值,mu_s是太阳方向天顶角的cos值,
 *  nu是向量pq与太阳单位方向向量的夹角cos值,ray_r_mu_intersects_ground射线是否与地面相交
 **/
void ComputeSingleScattering(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground,
	OUT(IrradianceSpectrum)rayleigh, OUT(IrradianceSpectrum)mie ){
	/* 数值积分采样个数 */
	const int SAMPLE_COUNT = 50;
	/* 积分步长,取一个最近的交点作为积分终点,或地面或大气顶层 */
	Length dx =
		DistanceToNearestAtmosphereBoundary(atmosphere, r, mu,
						     ray_r_mu_intersects_ground)/Number( SAMPLE_COUNT);
	DimensionlessSpectrum	rayleigh_sum	= DimensionlessSpectrum(0.0);
	DimensionlessSpectrum	mie_sum		= DimensionlessSpectrum(0.0);
	for ( int i = 0; i <= SAMPLE_COUNT; ++i ){
		Length d_i = Number(i)*dx;
		/* 当前点的rayleigh散射系数与mie散射系数. */
		DimensionlessSpectrum	rayleigh_i;
		DimensionlessSpectrum	mie_i;
		ComputeSingleScatteringIntegrand(atmosphere, transmittance_texture,
						  r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
		Number weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
		rayleigh_sum	  += rayleigh_i * weight_i;
		mie_sum		  += mie_i * weight_i;
	}
	// 此时尚未乘上相位函数,为了减少计算量,因为是均匀采样,所以将dx放到外循环去
	rayleigh = rayleigh_sum * dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering;
	mie 	  = mie_sum * dx * atmosphere.solar_irradiance * atmosphere.mie_scattering;
}

/**
 * 功能:
 *  计算rayleigh相位函数
 * 传入参数：
 *  nu是向量pq与太阳单位方向向量的夹角cos值
 **/
InverseSolidAngle RayleighPhaseFunction(Number nu){
	InverseSolidAngle k = 3.0 / (16.0 * PI * sr);// sr为立体角单位
	return(k * (1.0 + nu * nu) );
}


/**
 * 功能:
 *  计算mie相位函数
 * 传入参数：
 *  nu是向量pq与太阳单位方向向量的夹角cos值,g是散射的对称性因子
 *  g为正数表示光线大多数向后方散射,为负数说明更多的光线向前方散射
 **/
InverseSolidAngle MiePhaseFunction( Number g, Number nu ){
	InverseSolidAngle k = 3.0 / (8.0 * PI * sr) * (1.0 - g * g) / (2.0 + g * g);
	return(k * (1.0 + nu * nu) / pow( 1.0 + g * g - 2.0 * g * nu, 1.5 ) );
}

/**
 * 功能:
 *  将计算散射积分需要的四个参数(r,mu,mu_s,nu)映射到纹理坐标(u,v,w,z)
 * 传入参数：
 *  atmosphere大气参数模型,r为视点p高度,mu为视线天顶角cos值,
 *  mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值,
 *  ray_r_mu_intersects_ground射线是否与地面相交
 **/
vec4 GetScatteringTextureUvwzFromRMuMuSNu( IN(AtmosphereParameters)atmosphere,
					   Length r, Number mu, Number mu_s, Number nu,
					   bool ray_r_mu_intersects_ground ){
	/* 过视点的与地表相切的射线的切点到大气顶层的距离 */
	Length H = sqrt( atmosphere.top_radius * atmosphere.top_radius -
			 			atmosphere.bottom_radius * atmosphere.bottom_radius );
	/* 视点p到过视线的与地表相切的切线的切点的距离 */
	Length rho =
			SafeSqrt( r * r - atmosphere.bottom_radius * atmosphere.bottom_radius );
	Number u_r = GetTextureCoordFromUnitRange(rho/H, SCATTERING_TEXTURE_R_SIZE);

	/* 二元一次方程判别式,用于二元一次方程求根公式(求解射线(r,mu)与地表的交点) */
	Length	r_mu = r * mu;
	Area	discriminant	=
				r_mu*r_mu - r*r+atmosphere.bottom_radius*atmosphere.bottom_radius;
	Number u_mu;
	/* 若射线(r,mu)与地面有交点 */
	if ( ray_r_mu_intersects_ground ){
		/* 射线(r,mu)射点到与地面交点的距离d,纯粹的求根公式,以及d的上、下界 */
		Length d 		= -r_mu - SafeSqrt(discriminant);
		Length	d_min	= r - atmosphere.bottom_radius;
		Length	d_max	= rho;
		u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(
								 d_max == d_min ? 0.0:(d - d_min) / (d_max - d_min),
								 SCATTERING_TEXTURE_MU_SIZE / 2 );
	} else {
		/* 计算射线(r,mu)射点到大气层顶层边界交点的距离及其上界、下界 */
		Length	d		= -r_mu + SafeSqrt( discriminant + H * H );
		Length	d_min	= atmosphere.top_radius - r;
		Length	d_max	= rho + H;
		u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
								(d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE/2);
	}
	
	/* 对于mu_s,借用地表上有相同mu_s的点进行映射 */
	Length d 		  = DistanceToTopAtmosphereBoundary(
								atmosphere, atmosphere.bottom_radius, mu_s );
	Length	d_min	  = atmosphere.top_radius-atmosphere.bottom_radius;
	Length	d_max	  = H;
	Number	a		  = (d - d_min) / (d_max - d_min);
	Number	A		  = -2.0*atmosphere.mu_s_min*atmosphere.bottom_radius/(d_max-d_min);
	Number	u_mu_s   = GetTextureCoordFromUnitRange(max(1.0-a/A, 0.0)/(1.0+a),
							  SCATTERING_TEXTURE_MU_S_SIZE);
	Number u_nu = (nu + 1.0) / 2.0; //将nu从[-1,1]映射到[0,1]
	return vec4(u_nu, u_mu_s, u_mu, u_r);
}

/**
 * 功能:
 *  将纹理坐标(u,v,w,z)映射到计算散射积分需要的四个参数(r,mu,mu_s,nu)
 * 传入参数：
 *  atmosphere大气参数模型,r为视点p高度,mu为视线天顶角cos值,uvwz是4D纹理坐标
 *  mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值,
 *  ray_r_mu_intersects_ground射线是否与地面相交
 **/
void GetRMuMuSNuFromScatteringTextureUvwz( IN(AtmosphereParameters)atmosphere,
					   IN(vec4)uvwz, OUT(Length)r, OUT(Number)mu, OUT(Number)mu_s,
					   OUT(Number)nu, OUT(bool)ray_r_mu_intersects_ground ){
	Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius -
			 		   atmosphere.bottom_radius * atmosphere.bottom_radius );
	Length rho = H * GetUnitRangeFromTextureCoord( uvwz.w, SCATTERING_TEXTURE_R_SIZE );
	r = sqrt( rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius );
	if (uvwz.z < 0.5){
		Length	d_min	= r - atmosphere.bottom_radius;
		Length	d_max	= rho;
		Length	d	= d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
						1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2 );
		mu = d == 0.0 * m ? Number( -1.0 ) :
		     			ClampCosine( -(rho * rho + d * d) / (2.0 * r * d) );
		ray_r_mu_intersects_ground = true;
	} else {
		Length	d_min	= atmosphere.top_radius - r;
		Length	d_max	= rho + H;
		Length	d	= d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
						2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2 );
		mu = d == 0.0 * m ? Number( 1.0 ) :
		     			ClampCosine( (H * H - rho * rho - d * d) / (2.0 * r * d) );
		ray_r_mu_intersects_ground = false;
	}
	Number x_mu_s =
						GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
	Length	d_min	= atmosphere.top_radius - atmosphere.bottom_radius;
	Length	d_max	= H;
	Number	A	= -2.0 * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
	Number	a	= (A - x_mu_s * A) / (1.0 + x_mu_s * A);
	Length	d	= d_min + min( a, A ) * (d_max - d_min);
	mu_s = d == 0.0 * m ? Number( 1.0 ) :
	       			ClampCosine( (H * H - d * d) / (2.0 * atmosphere.bottom_radius * d) );
	nu = ClampCosine( uvwz.x * 2.0 - 1.0 );
}

/**
 * 功能:
 *  实际中只有3D纹理坐标,这里将3D转成4D,然后根据4D纹理获取(r,mu,mu_s,nu)参数
 * 传入参数：
 *  atmosphere大气参数模型,r为视点p高度,mu为视线天顶角cos值,uvwz是4D纹理坐标
 *  mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值,
 *  ray_r_mu_intersects_ground射线是否与地面相交
 **/
void GetRMuMuSNuFromScatteringTextureFragCoord(
	IN(AtmosphereParameters)atmosphere, IN(vec3)frag_coord,
	OUT(Length)r, OUT(Number)mu, OUT(Number)mu_s, OUT(Number)nu,
	OUT(bool)ray_r_mu_intersects_ground ){
	const vec4 SCATTERING_TEXTURE_SIZE = vec4(
						SCATTERING_TEXTURE_NU_SIZE - 1,
						SCATTERING_TEXTURE_MU_S_SIZE,
						SCATTERING_TEXTURE_MU_SIZE,
						SCATTERING_TEXTURE_R_SIZE);
	/* nu和mu_s的坐标值均由frag_cood.x获取,前者取整,后者取模 */
	Number frag_coord_nu =
				floor(frag_coord.x/Number(SCATTERING_TEXTURE_MU_S_SIZE));
	Number frag_coord_mu_s =
				mod(frag_coord.x, Number(SCATTERING_TEXTURE_MU_S_SIZE));
	vec4 uvwz = vec4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z)/
						SCATTERING_TEXTURE_SIZE;
	/* 根据uvwz这个4D纹理坐标从散射纹理中获取对应的(r,mu,mu_s,nu) */
	GetRMuMuSNuFromScatteringTextureUvwz(
				atmosphere, uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground );
	/* 对于nu,根据给定的mu和mu_s对其做一些上下界约束[cos(x+y),cos(x-y)] */
	nu = clamp(nu, mu * mu_s - sqrt( (1.0 - mu * mu) * (1.0 - mu_s * mu_s) ),
		    			mu * mu_s + sqrt( (1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

/**
 * 功能:
 *  有了以上的函数,现在我们可以计算一个指定的纹理单元对应的单次散射,结果存入rayleigh和mie
 * 传入参数：
 *  atmosphere大气参数模型,frag_coord为3D纹理坐标
 *  ray_r_mu_intersects_ground射线是否与地面相交
 **/
void ComputeSingleScatteringTexture(IN(AtmosphereParameters)atmosphere,
				     IN(TransmittanceTexture)transmittance_texture, IN(vec3)frag_coord,
				     OUT(IrradianceSpectrum)rayleigh, OUT(IrradianceSpectrum)mie){
	Length	r;
	Number	mu;
	Number	mu_s;
	Number	nu;
	bool	ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
						   r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	ComputeSingleScattering(atmosphere, transmittance_texture,
				 r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

/**
 * 功能:
 *  获取射点到与最近大气层边界交点之间的散射量,需要两次3D的纹理查找
 * 传入参数：
 *  atmosphere大气参数模型,scattering_texture为散射预存纹理,r为视点p海拔高度,
 *  mu为视线天顶角cos值,mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值
 *  ray_r_mu_intersects_ground射线是否与地面相交
 **/
AbstractSpectrum GetScattering(IN(AtmosphereParameters)atmosphere,
	IN(AbstractScatteringTexture)scattering_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground ){
	/* 根据给定的(r,mu,mu_s,nu)计算对应的4D纹理坐标uvwz */
	vec4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
						atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground );
	Number	tex_coord_x	= uvwz.x * Number(SCATTERING_TEXTURE_NU_SIZE - 1);
	Number	tex_x			= floor(tex_coord_x); /* 整数部分 */
	Number	lerp			= tex_coord_x-tex_x;  /* 小数部分 */
	vec3 uvw0 = vec3((tex_x + uvwz.y) / Number( SCATTERING_TEXTURE_NU_SIZE ),
			  			uvwz.z, uvwz.w);
	vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / Number( SCATTERING_TEXTURE_NU_SIZE ),
			  			uvwz.z, uvwz.w);
	/* 根据lerp线性插值 */
	return AbstractSpectrum(texture(scattering_texture, uvw0) * (1.0 - lerp) +
				 			    texture(scattering_texture, uvw1) * lerp);
}

/**
 * 功能:
 *  获取射点到与大气层边界交点(或大气顶层,或地表面)之间的散射辐照度,需要两次3D的纹理查找
 * 传入参数：
 *  atmosphere大气参数模型,single_rayleigh_scattering_texture为rayleigh单次散射纹理,
 *  single_mie_scattering_texture为mie单词散射纹理,multiple_scattering_texture
 *  r为视点p高度,mu为视线天顶角cos值,mu_s是太阳方向天顶角的cos值,
 *  nu是向量pq与太阳单位方向向量的夹角cos值,ray_r_mu_intersects_ground射线是否与地面相交
 *  multiple_scattering_texture为多次散射纹理,scattering_order为散射重数
 **/
RadianceSpectrum GetScattering(
    IN(AtmosphereParameters) atmosphere,
    IN(ReducedScatteringTexture) single_rayleigh_scattering_texture,
    IN(ReducedScatteringTexture) single_mie_scattering_texture,
    IN(ScatteringTexture) multiple_scattering_texture,
    Length r, Number mu, Number mu_s, Number nu,
    bool ray_r_mu_intersects_ground,
    int scattering_order) {
	if (scattering_order == 1) {//单次散射
		IrradianceSpectrum rayleigh = GetScattering(
        		atmosphere, single_rayleigh_scattering_texture, r, mu, mu_s, nu,
        		ray_r_mu_intersects_ground);
		IrradianceSpectrum mie = GetScattering(
        		atmosphere, single_mie_scattering_texture, r, mu, mu_s, nu,
        		ray_r_mu_intersects_ground);
	 // 最后还要乘上相应的相位函数
    return rayleigh * RayleighPhaseFunction(nu) +
        		mie * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
	} else {//多次散射
	return GetScattering(atmosphere, multiple_scattering_texture, r, mu, mu_s, nu,
        						ray_r_mu_intersects_ground);
  }
}

/**
 * 功能:
 *  经过n-2次散射地面接收的辐照度
 * 传入参数：
 *  atmosphere大气参数模型,irradiance_texture辐照度纹理,
 *  single_mie_scattering_texture为mie单词散射纹理,multiple_scattering_texture
 *  r为视点p海拔高度,mu为视线天顶角cos值,mu_s是太阳方向天顶角的cos值,
 **/
IrradianceSpectrum GetIrradiance(IN(AtmosphereParameters)atmosphere,
									   IN(IrradianceTexture)irradiance_texture,
										Length r, Number mu_s);

/**
 * 功能:
 *  根据n-1重的散射纹理计算第n重的散射密度,返回辐射密度谱
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为透射率纹理,
 *  single_rayleigh_scattering_texture单次rayleigh散射纹理,
 *  single_mie_scattering_texture为mie单次散射纹理,multiple_scattering_texture多次散射纹理
 *  irradiance_texture地面接收的辐照度纹理,r为积分点p高度,mu为视线天顶角cos值,
 *  mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值,
 *  scattering_order散射次数
 **/
RadianceDensitySpectrum ComputeScatteringDensity(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ReducedScatteringTexture)single_rayleigh_scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	IN(ScatteringTexture)multiple_scattering_texture,
	IN(IrradianceTexture)irradiance_texture,
	Length r, Number mu, Number mu_s, Number nu, int scattering_order ){
	vec3 zenith_direction = vec3(0.0, 0.0, 1.0);
	vec3	omega		= vec3(sqrt(1.0 - mu * mu), 0.0, mu);/* 视线方向向量 */
	Number	sun_dir_x	= omega.x == 0.0 ? 0.0 : (nu - mu * mu_s) / omega.x;
	Number	sun_dir_y	= sqrt( max( 1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0 ) );
	vec3	omega_s	= vec3(sun_dir_x, sun_dir_y, mu_s);/* 太阳方向向量 */
	const int		SAMPLE_COUNT	= 16;//积分采样数u
	const Angle		dphi		= pi / Number( SAMPLE_COUNT );//采样步长
	const Angle		dtheta		= pi / Number( SAMPLE_COUNT );
	RadianceDensitySpectrum rayleigh_mie	=
			RadianceDensitySpectrum( 0.0 * watt_per_cubic_meter_per_sr_per_nm );
	/* 双重积分,theta为天顶角 */
	for (int l = 0; l < SAMPLE_COUNT; ++l){
		Angle	theta		= (Number( l ) + 0.5) * dtheta;
		Number	cos_theta	= cos( theta );
		Number	sin_theta	= sin( theta );
		/* 判断是否与地面有交点 */
		bool ray_r_theta_intersects_ground = RayIntersectsGround(atmosphere, r, cos_theta);
		/* 当前射线到地表交点的距离以及光学深度仅取决于theta(即天顶角),所以放到外循环计算 */
		Length	distance_to_ground = 0.0 * m;
		DimensionlessSpectrum transmittance_to_ground = DimensionlessSpectrum( 0.0 );
		DimensionlessSpectrum ground_albedo = DimensionlessSpectrum( 0.0 );
		/* 与地面有交点 */
		if (ray_r_theta_intersects_ground){
			/* 射点沿射线到地面的距离 */
			distance_to_ground = DistanceToBottomAtmosphereBoundary(atmosphere, r, cos_theta);
			/* 相应的到地面的透射率 */
			transmittance_to_ground =
				GetTransmittance( atmosphere, transmittance_texture, r, cos_theta,
						  distance_to_ground, true /* ray_intersects_ground */ );
			ground_albedo = atmosphere.ground_albedo;/* 地面反照率 */
		}
		for ( int m = 0; m < 2 * SAMPLE_COUNT; ++m ){/* 内循环积分,这里积的是方位角 */
			Angle phi = (Number(m) + 0.5) * dphi;
			/* 由phi(方位角)和theta(天顶角)两个角度指定的方向向量 */
			vec3 omega_i = vec3( cos( phi ) * sin_theta, sin( phi ) * sin_theta, cos_theta );
			/* 立体角domega=sin(theta)*dtheta*dphi; */
			SolidAngle domega_i = (dtheta / rad) * (dphi / rad) * sin( theta ) * sr;
			/*经过n-1次反射从omega_i方向的辐射度L_i由n-1重散射的散射值累加。
			 * 太阳方向向量与omge_i向量夹角的cos值 */
			Number nu1 = dot(omega_s, omega_i);
			/* 获取n-1散射在此采样点接收的入射光线辐射度 */
			RadianceSpectrum incident_radiance = GetScattering(atmosphere,
									    single_rayleigh_scattering_texture, single_mie_scattering_texture,
									    multiple_scattering_texture, r, omega_i.z, mu_s, nu1,
									    ray_r_theta_intersects_ground, scattering_order - 1);
			/* 采样点接收地面辐射的光线,这部分的值主要是由采样点到地面的光学深度、
			 * 地面反照率、地面的BRDF以及地面吸收的n-2次反射的辐照度的乘积 */
			vec3 ground_normal = normalize(zenith_direction * r + omega_i * distance_to_ground);
			/* 地面吸收的n-2次反射的辐照度 */
			IrradianceSpectrum ground_irradiance = GetIrradiance(
								atmosphere, irradiance_texture, atmosphere.bottom_radius,
								dot(ground_normal, omega_s));
			incident_radiance += transmittance_to_ground *
					     ground_albedo * (1.0 / (PI * sr)) * ground_irradiance;

			/* 从omega_i方向向-omega方向散射的辐射度为incident_radiance、
			 * 散射系数、方向omega和omega_i的相位函数之积 */
			Number	nu2			= dot(omega, omega_i);
			Number	rayleigh_density	= GetProfileDensity(
						atmosphere.rayleigh_density, r - atmosphere.bottom_radius );
			Number mie_density = GetProfileDensity(
						atmosphere.mie_density, r - atmosphere.bottom_radius );
			rayleigh_mie += incident_radiance * (
						atmosphere.rayleigh_scattering * rayleigh_density *
						RayleighPhaseFunction( nu2 ) +
						atmosphere.mie_scattering * mie_density *
						MiePhaseFunction( atmosphere.mie_phase_function_g, nu2 ) ) * domega_i;
		}
	}
	return rayleigh_mie;
}

/**
 * 功能:
 *  计算多重散射系数(从视点r处到与边界(或大气顶层或地面)交点的内散射积分)
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为透射率纹理,
 *  scattering_density_texture散射密度纹理,
 *  r为视点p高度,mu为视线天顶角cos值,
 *  mu_s是太阳方向天顶角的cos值,nu是向量pq与太阳单位方向向量的夹角cos值,
 *  ray_r_mu_intersects_ground是否与地面相交
 **/
RadianceSpectrum ComputeMultipleScattering(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ScatteringDensityTexture)scattering_density_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground ){
	const int SAMPLE_COUNT = 50; /* 积分采样数 */
	Length dx =	/* 积分步长 */
				DistanceToNearestAtmosphereBoundary(
					atmosphere, r, mu, ray_r_mu_intersects_ground )/Number(SAMPLE_COUNT);
	RadianceSpectrum rayleigh_mie_sum =
				RadianceSpectrum( 0.0 * watt_per_square_meter_per_sr_per_nm );
	for ( int i = 0; i <= SAMPLE_COUNT; ++i ){
		Length d_i = Number( i ) * dx;
		/* 当前采样积分点的高度r_i */
		Length r_i = ClampRadius(atmosphere, sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r));
		Number	mu_i	= ClampCosine((r * mu + d_i) / r_i);
		Number	mu_s_i	= ClampCosine((r * mu_s + d_i * nu) / r_i);
		/* 当前采样点的rayleigh散射系数、mie散射系数 */
		RadianceSpectrum rayleigh_mie_i =
				GetScattering(
					atmosphere, scattering_density_texture, r_i, mu_i, mu_s_i, nu,
					ray_r_mu_intersects_ground ) *
				GetTransmittance(
					atmosphere, transmittance_texture, r, mu, d_i,
					ray_r_mu_intersects_ground ) * dx;
		Number weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
		rayleigh_mie_sum += rayleigh_mie_i * weight_i;
	}
	return rayleigh_mie_sum;
}

/**
 * 功能:
 *  计算当前3D纹理坐标对应的(r,mu,mu_s,nu)处接收的散射辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,
 *  single_rayleigh_scattering_texture单次rayleigh散射纹理,
 *  single_mie_scattering_texture为mie单次散射纹理,
 *  multiple_scattering_texture多次散射纹理
 *  irradiance_texture地面接收的辐照度纹理,frag_coord为3D纹理坐标,
 *  scattering_order散射重数
 **/
RadianceDensitySpectrum ComputeScatteringDensityTexture(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ReducedScatteringTexture)single_rayleigh_scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	IN(ScatteringTexture)multiple_scattering_texture,
	IN(IrradianceTexture)irradiance_texture,
	IN(vec3)frag_coord, int scattering_order ){
	Length	r;
	Number	mu;
	Number	mu_s;
	Number	nu;
	bool	ray_r_mu_intersects_ground;
	/* 从frag_coord的3D纹理坐标计算出当前的(r,mu,mu_s,nu) */
	GetRMuMuSNuFromScatteringTextureFragCoord( atmosphere, frag_coord,
						   r, mu, mu_s, nu, ray_r_mu_intersects_ground );
	return ComputeScatteringDensity(atmosphere, transmittance_texture,
					 single_rayleigh_scattering_texture, single_mie_scattering_texture,
					 multiple_scattering_texture, irradiance_texture, r, mu, mu_s, nu,
					 scattering_order);
}


/**
 * 功能:
 *  计算当前3D纹理坐标对应的(r,mu,mu_s,nu)处的多重散射系数
 *  从视点到大气顶层的路径上的光学深度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,
 *  scattering_density_texture为散射谱纹理,frag_coord为3D纹理坐标,
 *  nu为射线(r,mu)与太阳方向向量夹角cos值
 **/
RadianceSpectrum ComputeMultipleScatteringTexture(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ScatteringDensityTexture)scattering_density_texture,
	IN(vec3)frag_coord, OUT(Number)nu ){
	Length	r;
	Number	mu;
	Number	mu_s;
	bool	ray_r_mu_intersects_ground;
	/* 从frag_coord的3D纹理坐标计算出当前的(r,mu,mu_s,nu) */
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
						   r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	return ComputeMultipleScattering(atmosphere, transmittance_texture,
					  scattering_density_texture, r, mu, mu_s, nu,
					  ray_r_mu_intersects_ground);
}

/**
 * 功能:
 *  计算地面接收的直接辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,
 *  scattering_density_texture为散射谱纹理,(r,mu_s)为(视点高度,天顶角)
 **/
IrradianceSpectrum ComputeDirectIrradiance(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	Length r, Number mu_s){
	Number alpha_s = atmosphere.sun_angular_radius / rad;//太阳角半径
	/* 太阳圆盘可见部分的近似平均余弦因子 */
	Number average_cosine_factor = mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
					 (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s) );
	return atmosphere.solar_irradiance *
	       GetTransmittanceToTopAtmosphereBoundary(atmosphere,
				 transmittance_texture, r, mu_s ) * average_cosine_factor;
}

/**
 * 功能:
 *  计算地面的间接辐照度,对以地面法线为轴的半球进行积分
 * 传入参数：
 *  atmosphere大气参数模型,single_rayleigh_scattering_texture为单次rayleigh散射纹理,
 *  single_mie_scattering_texture为单次mie散射纹理,
 *  multiple_scattering_texture多次散射纹理,
 *  (r,mu_s)为(海拔高度,天顶角),scattering_order散射次数
 **/
IrradianceSpectrum ComputeIndirectIrradiance(
	IN(AtmosphereParameters)atmosphere,
	IN(ReducedScatteringTexture)single_rayleigh_scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	IN(ScatteringTexture)multiple_scattering_texture,
	Length r, Number mu_s, int scattering_order ){
	const int	SAMPLE_COUNT	= 32;/* 积分采样数 */
	const Angle	dphi		= pi / Number( SAMPLE_COUNT );/* 采样步长 */
	const Angle	dtheta		= pi / Number( SAMPLE_COUNT );

	IrradianceSpectrum result =
			IrradianceSpectrum( 0.0 * watt_per_square_meter_per_nm );
	/* 射线(r,mu_s)的方向向量 */
	vec3 omega_s = vec3( sqrt( 1.0 - mu_s * mu_s ), 0.0, mu_s );
	/* 对整个半球进行双重积分 */
	for ( int j = 0; j < SAMPLE_COUNT / 2; ++j ){//这里除以2是因为只对半球积分
		Angle theta = (Number( j ) + 0.5) * dtheta;
		for ( int i = 0; i < 2 * SAMPLE_COUNT; ++i ){
			Angle phi = (Number( i ) + 0.5) * dphi;
			vec3 omega = /* (theta,phi)指定的方向向量 */
				vec3( cos( phi ) * sin( theta ), sin( phi ) * sin( theta ), cos( theta ) );
			/* 立体角微元domega */
			SolidAngle domega = (dtheta / rad) * (dphi / rad) * sin( theta ) * sr;
			/* omega与omega_s的夹角cos值 */
			Number nu = dot( omega, omega_s );
			result += GetScattering( atmosphere, single_rayleigh_scattering_texture,
						 single_mie_scattering_texture, multiple_scattering_texture,
						 r, omega.z, mu_s, nu, false /* ray_r_theta_intersects_ground */,
						 scattering_order ) * omega.z * domega;
		}
	}
	return result;
}

/**
 * 功能:
 *  将参数(r,mu_s)映射到(u,v)纹理坐标,用于地面接收的辐照度的预计算
 *  因为地面辐照度的计算仅考虑了水平面,所以辐照度计算仅需参数(r,mu_s)
 * 传入参数：
 *  atmosphere大气参数模型, (r,mu_s)为(高度,天顶角)
 **/
vec2 GetIrradianceTextureUvFromRMuS( IN(AtmosphereParameters)atmosphere,
				     Length r, Number mu_s ){
	Number x_r = (r - atmosphere.bottom_radius) /
		     (atmosphere.top_radius - atmosphere.bottom_radius);/* 将r投影到[0,1]之间 */
	Number x_mu_s = mu_s * 0.5 + 0.5;/* 同理将mu_s由[-1,1]投影到[0,1] */
	return(vec2( GetTextureCoordFromUnitRange( x_mu_s, IRRADIANCE_TEXTURE_WIDTH ),
		     GetTextureCoordFromUnitRange( x_r, IRRADIANCE_TEXTURE_HEIGHT ) ) );
}

/**
 * 功能:
 *  将纹理坐标(u,v)映射到参数(r,mu_s),GetIrradianceTextureUvFromRMuS的逆过程
 * 传入参数：
 *  atmosphere大气参数模型,uv为2D纹理坐标,(r,mu_s)为(高度,天顶角)
 **/
void GetRMuSFromIrradianceTextureUv( IN(AtmosphereParameters)atmosphere,
				     IN(vec2)uv, OUT(Length)r, OUT(Number)mu_s ){
	Number	x_mu_s	= GetUnitRangeFromTextureCoord( uv.x, IRRADIANCE_TEXTURE_WIDTH );
	Number	x_r		= GetUnitRangeFromTextureCoord( uv.y, IRRADIANCE_TEXTURE_HEIGHT );
	r = atmosphere.bottom_radius +
	    	x_r * (atmosphere.top_radius - atmosphere.bottom_radius);
	mu_s = ClampCosine( 2.0 * x_mu_s - 1.0 );
}

// 辐照度纹理大小
const vec2 IRRADIANCE_TEXTURE_SIZE =
    vec2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

/**
 * 功能:
 *  计算地面的直接从太阳接收的辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,
 *  frag_coord为2D纹理坐标
 **/
IrradianceSpectrum ComputeDirectIrradianceTexture(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(vec2)frag_coord ){
	Length	r;
	Number	mu_s;
	/* 根据纹理坐标frag_coord计算出(r,mu_s) */
	GetRMuSFromIrradianceTextureUv(
		atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s );
	/* 用于计算直接辐照度 */
	return(ComputeDirectIrradiance( atmosphere, transmittance_texture, r, mu_s ) );
}

/**
 * 功能:
 *  计算地面的间接辐照度
 * 传入参数：
 *  atmosphere大气参数模型,single_rayleigh_scattering_texture单重rayleigh散射纹理,
 *  single_mie_scattering_texture为单重mie散射纹理,
 *  multiple_scattering_texture为多重散射纹理
 *  frag_coord为2D纹理坐标,scattering_order为散射重数
 **/
IrradianceSpectrum ComputeIndirectIrradianceTexture(
	IN(AtmosphereParameters)atmosphere,
	IN(ReducedScatteringTexture)single_rayleigh_scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	IN(ScatteringTexture)multiple_scattering_texture,
	IN(vec2)frag_coord, int scattering_order ){
	Length	r;
	Number	mu_s;
	GetRMuSFromIrradianceTextureUv(
			atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s );
	return(ComputeIndirectIrradiance( atmosphere,single_rayleigh_scattering_texture,
					 	single_mie_scattering_texture,multiple_scattering_texture,
						 r, mu_s, scattering_order ) );
}

/**
 * 功能:
 *  通过对辐照度纹理查找一次,获取地面辐照度值
 * 传入参数：
 *  atmosphere大气参数模型,irradiance_texture为辐照度纹理,
 *  (r,mu_s)为(高度,天顶角)
 **/
IrradianceSpectrum GetIrradiance(
	IN(AtmosphereParameters)atmosphere,
	IN(IrradianceTexture)irradiance_texture,
	Length r, Number mu_s ){
	vec2 uv = GetIrradianceTextureUvFromRMuS(atmosphere, r, mu_s);
	return IrradianceSpectrum(texture( irradiance_texture, uv));
}

/**
 * 功能:
 *  在将rayleigh和mie散射纹理合并时,从中推出mie单重散射值
 * 传入参数：
 *  atmosphere大气参数模型,scattering为合并的单重散射值
 **/
#ifdef COMBINED_SCATTERING_TEXTURES
vec3 GetExtrapolatedSingleMieScattering(
	IN(AtmosphereParameters)atmosphere, IN(vec4)scattering ){
	if (scattering.r == 0.0)return vec3(0.0);
	return scattering.rgb * scattering.a/scattering.r *
	      (atmosphere.rayleigh_scattering.r/atmosphere.mie_scattering.r) *
	      (atmosphere.mie_scattering / atmosphere.rayleigh_scattering);
}
#endif

/**
 * 功能:
 *  从纹理中获取内散射值(在合并纹理的情况调用以下函数)
 * 传入参数：
 *  atmosphere大气参数模型,scattering_texture散射纹理,
 *  single_mie_scattering_texture为单次mie散射纹理,
 *  (r,mu,mu_s,nu)为(高度,天顶角,太阳天顶角,射线(r,mu)与太阳方向向量夹角cos值)
 *  ray_r_mu_intersects_ground是否与地面相交,
 *  single_mie_scattering为单重mie散射辐照度
 **/
IrradianceSpectrum GetCombinedScattering(
	IN(AtmosphereParameters)atmosphere,
	IN(ReducedScatteringTexture)scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	Length r, Number mu, Number mu_s, Number nu,
	bool ray_r_mu_intersects_ground,
	OUT(IrradianceSpectrum)single_mie_scattering){
	/* 将(nu,mu_s,mu,r)映射到纹理坐标(u,v,w,z),与前面一样,还要将4D纹理坐标映射到3D */
	vec4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
			atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground );
	Number	tex_coord_x	= uvwz.x * Number(SCATTERING_TEXTURE_NU_SIZE - 1);
	Number	tex_x			= floor( tex_coord_x ); /* 整数部分 */
	Number	lerp			= tex_coord_x - tex_x;  /* 小数部分 */
	vec3	uvw0			= vec3( (tex_x + uvwz.y) / Number( SCATTERING_TEXTURE_NU_SIZE ),
									uvwz.z, uvwz.w );
	vec3	uvw1 			= vec3( (tex_x + 1.0 + uvwz.y) / Number( SCATTERING_TEXTURE_NU_SIZE ),
			  						uvwz.z, uvwz.w );
	/* 对于合并成一个纹理的方案,只采样scattering_texture,然后分出其中的mie_scattering */
#ifdef COMBINED_SCATTERING_TEXTURES
	vec4 combined_scattering =
			texture( scattering_texture, uvw0 ) * (1.0 - lerp) +
			texture( scattering_texture, uvw1 ) * lerp;
	IrradianceSpectrum scattering = IrradianceSpectrum( combined_scattering );
	single_mie_scattering =
			GetExtrapolatedSingleMieScattering( atmosphere, combined_scattering );
#else/* 非合并的情况,直接采样各自的纹理 */
	IrradianceSpectrum scattering = IrradianceSpectrum(
			texture( scattering_texture, uvw0 ) * (1.0 - lerp) +
			texture( scattering_texture, uvw1 ) * lerp );
	single_mie_scattering = IrradianceSpectrum(
			texture( single_mie_scattering_texture, uvw0 ) * (1.0 - lerp) +
			texture( single_mie_scattering_texture, uvw1 ) * lerp );
#endif
	return scattering;
}

/**
 * 功能:
 *  获取天空的辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,scattering_texture散射纹理,
 *  single_mie_scattering_texture为单次mie散射纹理,camera视点位置,view_ray视线向量
 *  shadow_length阴影长度,由阴影体算法计算得到,sun_direction太阳方向向量,transmittance透射率
 **/
RadianceSpectrum GetSkyRadiance(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ReducedScatteringTexture)scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	Position camera, IN(Direction)view_ray, Length shadow_length,
	IN(Direction)sun_direction, OUT(DimensionlessSpectrum)transmittance){
	Length r = length(camera);/* 视点所在的高度 */
	Length rmu = dot( camera, view_ray );/* r*天顶角cos值 */
	/* 视点沿视线到大气层顶层的距离 */
	Length distance_to_top_atmosphere_boundary = -rmu -
					sqrt( rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius );
	/* 如果观察者在太空且视线与大气层有交点,把观察者移到视线与大气层顶部交点的位置 */
	if ( distance_to_top_atmosphere_boundary > 0.0 * m ){
		camera	= camera + view_ray * distance_to_top_atmosphere_boundary;
		r	= atmosphere.top_radius;
		rmu	+= distance_to_top_atmosphere_boundary;
	} else if (r > atmosphere.top_radius){
		/* 视线与大气层无交点,直接返回0 */
		transmittance = DimensionlessSpectrum( 1.0 );
		return(RadianceSpectrum( 0.0 * watt_per_square_meter_per_sr_per_nm));
	}
	/* 计算(r,mu,mu_s,nu)即(高度,天顶角cos,太阳天顶角cos,太阳方向向量与视线夹角cos) */
	Number	mu		= rmu / r;
	Number	mu_s	= dot( camera, sun_direction ) / r;
	Number	nu		= dot( view_ray, sun_direction );
	/* 检测视线是否与地面相交 */
	bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);
	/* 与地面相交则光学深度为0 */
	transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum( 0.0 ) :
			GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r, mu );
	IrradianceSpectrum	single_mie_scattering;
	IrradianceSpectrum	scattering;
	if ( shadow_length == 0.0 * m ){/* 如果不需要体积光效果 */
		scattering = GetCombinedScattering(
				atmosphere, scattering_texture, single_mie_scattering_texture,
				r, mu, mu_s, nu, ray_r_mu_intersects_ground,single_mie_scattering );
	} else {
		/* 实现体积光效果:我们省去从摄像机出发沿着视线长度为shadow_length的这段的散射计算, 
		 * 只计算剩下的那段,即沿视线长度为d处的点到大气顶层交点这一段 */
		Length	d		= shadow_length;
		Length	r_p		= ClampRadius(atmosphere, sqrt(d*d+2.0*r*mu*d+r*r));
		Number	mu_p	= (r * mu + d) / r_p;
		Number	mu_s_p	= (r * mu_s + d * nu) / r_p;
		scattering 	= GetCombinedScattering(
				atmosphere, scattering_texture, single_mie_scattering_texture,
				r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
				single_mie_scattering);
		/* 视点p到沿视线在长度shadow_length处的点的光学深度 */
		DimensionlessSpectrum shadow_transmittance =
					GetTransmittance(atmosphere, transmittance_texture,
					  		r, mu, shadow_length, ray_r_mu_intersects_ground);
		scattering = scattering * shadow_transmittance;
		single_mie_scattering = single_mie_scattering * shadow_transmittance;
	}
	return scattering * RayleighPhaseFunction(nu)
		 + single_mie_scattering * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

/**
 * 功能:
 *  获取视点到某一点的辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学深度纹理,scattering_texture散射纹理,
 *  single_mie_scattering_texture为单次mie散射纹理,camera视点位置,point目标位置,
 *  shadow_length阴影长度,sun_direction太阳方向向量,transmittance光学深度
 **/
RadianceSpectrum GetSkyRadianceToPoint(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(ReducedScatteringTexture)scattering_texture,
	IN(ReducedScatteringTexture)single_mie_scattering_texture,
	Position camera, IN(Position)point, Length shadow_length,
	IN(Direction)sun_direction, OUT(DimensionlessSpectrum)transmittance ){
	Direction	view_ray	= normalize( point - camera );/* 视线方向向量 */
	Length		r			= length(camera);             /* 海拔高度 */
	Length		rmu			= dot(camera, view_ray);      /* r*天顶角cos值 */
	/* 到达大气层顶部的距离 */
	Length distance_to_top_atmosphere_boundary = -rmu -
				sqrt(rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius);
	/* 如果视点在太空中, 且视线与大气层有交点,那么把view沿着视线移动到大气层顶部 */
	if ( distance_to_top_atmosphere_boundary > 0.0 * m ){
		camera	= camera + view_ray * distance_to_top_atmosphere_boundary;
		r		= atmosphere.top_radius;
		rmu	   += distance_to_top_atmosphere_boundary;
	}
	/* 计算(r,mu,mu_s,nu)参数用于第一次纹理查找,得到camera到大气层边界的内散射积分 */
	Number	mu		= rmu / r;
	Number	mu_s	= dot(camera, sun_direction)/r;
	Number	nu		= dot(view_ray, sun_direction);
	Length	d		= length(point - camera); /* 目标point到摄像机的距离 */
	/* 检测视线是否与地面相交 */
	bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);
	/* 获取在camera到point的光学长度 */
	transmittance = GetTransmittance(atmosphere, transmittance_texture,
					  		r, mu, d, ray_r_mu_intersects_ground);
	/* 获取camera到大气边界的rayleigh内散射积分 */
	IrradianceSpectrum	single_mie_scattering;
	IrradianceSpectrum	scattering = GetCombinedScattering(
			atmosphere, scattering_texture, single_mie_scattering_texture,
			r, mu, mu_s, nu, ray_r_mu_intersects_ground,
			single_mie_scattering);

	/* 计算(r,mu,mu_s,nu)用于第二次纹理查找,获取点point到大气层边界的内散射积分
	 * 如果需要实现体积光效果(shadow_length>0),那么我们应该忽略视线方向末端的shadow_length长度的散射
	 * 因此令d=d-shadow_length*/
	d = max(d - shadow_length, 0.0 * m);
	/* point处高度 */
	Length r_p = ClampRadius(atmosphere, sqrt( d * d + 2.0 * r * mu * d + r * r ));
	Number	mu_p	= (r * mu + d) / r_p;/* point处天顶角cos */
	Number	mu_s_p	= (r * mu_s + d * nu) / r_p;/* point处太阳方向天顶角cos */

	/* 查找point点处的内散射纹理 */
	IrradianceSpectrum	single_mie_scattering_p;
	IrradianceSpectrum	scattering_p = GetCombinedScattering(
			atmosphere, scattering_texture, single_mie_scattering_texture,
			r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
			single_mie_scattering_p );
	/* 把以上查找结果综合起来得到camera到point之间的散射 */
	DimensionlessSpectrum shadow_transmittance = transmittance;
	if ( shadow_length > 0.0 * m ){/* 要实现体积光效果 */
		shadow_transmittance = GetTransmittance( atmosphere, transmittance_texture,
							 		r, mu, d, ray_r_mu_intersects_ground );
	}
	/* camera到point的散射=camera到大气层边界的scattering-point到大气层边界的scattering */
	scattering = scattering - shadow_transmittance * scattering_p;
	single_mie_scattering	=
			single_mie_scattering - shadow_transmittance * single_mie_scattering_p;
#ifdef COMBINED_SCATTERING_TEXTURES
	/* 对于combined的方案,需要将mie单次散射从scattering提取出来 */
	single_mie_scattering = GetExtrapolatedSingleMieScattering(
				atmosphere, vec4(scattering, single_mie_scattering.r));
#endif
	/* 埃尔米特插值,避免太阳在水平线以下时失真 */
	single_mie_scattering = single_mie_scattering *
					smoothstep(Number(0.0), Number(0.01), mu_s);
	return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
	       MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

/**
 * 功能:
 *  计算地表的辐照度
 * 传入参数：
 *  atmosphere大气参数模型,transmittance_texture为光学长度纹理,irradiance_texture辐照度纹理,
 *  point目标位置,normal地面法线,sun_direction太阳方向向量,sky_irradiance天空s辐照度
 **/
IrradianceSpectrum GetSunAndSkyIrradiance(
	IN(AtmosphereParameters)atmosphere,
	IN(TransmittanceTexture)transmittance_texture,
	IN(IrradianceTexture)irradiance_texture,
	IN(Position)point, IN(Direction)normal, IN(Direction)sun_direction,
	OUT(IrradianceSpectrum)sky_irradiance ){
	Length r = length(point);/* point点海拔高度 */
	Number mu_s = dot(point, sun_direction) / r;/* 相应的天顶角 */
	/* 间接辐照度 */
	sky_irradiance = GetIrradiance( atmosphere, irradiance_texture, r, mu_s ) *
			 				(1.0 + dot( normal, point ) / r) * 0.5;
	/* 直接辐照度 */
	return atmosphere.solar_irradiance *
	       GetTransmittanceToSun(atmosphere, transmittance_texture, r, mu_s ) *
	       max(dot( normal, sun_direction), 0.0);
}