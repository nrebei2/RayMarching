#define PI 3.1415926535
// Sphere
// s: radius
float sdSphere(float3 p, float s) {
	return length(p) - s;
}

// Box
// b: size of box in x/y/z
float sdBox(float3 p, float3 b) {
	float3 d = abs(p) - b;
	return min(max(d.x, max(d.y, d.z)), 0.0) +
		length(max(d, 0.0));
}

// Rounded Box
float sdRoundBox(in float3 p, in float3 b, in float r) {
	float3 q = abs(p) - b;
	return min(max(q.x, max(q.y, q.z)), 0.0) + length(max(q, 0.0)) - r;
}

// Torus
float sdTorus(float3 p, float2 t) {
	float2 q = float2(length(p.xz)-t.x, p.y);
	return length(q) - t.y;
}

// (Infinite) Plane
// n.xyz: normal of the plane (normalized)
// n.w: offset
float sdPlane(float3 p, float4 n) {
	// n must be normalized
	return dot(p, n.xyz) + n.w;
}

// BOOLEAN OPERATORS //

// Union
float4 opU(float4 d1, float4 d2) {
	return (d1.w < d2.w) ? d1 : d2;
}

// Subtraction
float opS(float d1, float d2) {
	return max(-d1, d2);
}

// Intersection
float opI(float d1, float d2) {
	return max(d1, d2);
}

// Mod Position Axis
float pMod1 (inout float p, float size) {
	float halfsize = size * 0.5;
	float c = floor((p+halfsize)/size);
	p = fmod(p+halfsize,size)-halfsize;
	p = fmod(-p+halfsize,size)-halfsize;
	return c;
}

// SMOOTH BOOLEAN OPERATORS

// Union
float4 opUS(float4 d1, float4 d2, float k) {
	float h = clamp(0.5 + 0.5 * (d2.w - d1.w) / k, 0.0, 1.0);
	// Lerp between color and distance between objects
	float3 color = lerp(d2.rgb, d1.rgb, h);
	float dist =  lerp(d2.w, d1.w, h) - k * h * (1.0 - h);
	return float4(color, dist);
}

float4 printBulb(float4 d1, float k) {
	// Lerp between color and distance between objects
	float3 color = d1.rgb;
	float dist =  d1.w;
	return float4(color, dist);
}


// Subtraction
float opSS(float d1, float d2, float k) {
	float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
	return lerp(d2, -d1, h) + k * h * (1.0 - h);
}

// Intersection
float opIS(float d1, float d2, float k) {
	float h = clamp(0.5 - 0.5 * (d2 - d1) / k, 0.0, 1.0);
	return lerp(d2, d1, h) + k * h * (1.0 - h);
}

float4 qsqr( in float4 a ) // square a quaterion
{
    return float4( a.x*a.x - a.y*a.y - a.z*a.z - a.w*a.w,
                 2.0*a.x*a.y,
                 2.0*a.x*a.z,
                 2.0*a.x*a.w );
}

void sphereFold(inout float3 z, inout float dz)
{
	float r2 = dot(z,z);
	if (r2 < 0.5)
    { 
		float temp = 2.0;
		z *= temp;
		dz*= temp;
	}
    else if (r2 < 1.0)
    { 
		float temp = 1.0 / r2;
		z *= temp;
		dz*= temp;
	}
}

void boxFold(inout float3 z, inout float dz)
{
	z = clamp(z, -1.0, 1.0) * 2.0 - z;
}

//static mandelbulb
float sdMandelbulb(float3 p)
{
	float3 w = p;
    float m = dot(w, w);

	float dz = 1.0;
        
	for(int i = 0; i < 2; i++)
    {
        dz = 8 * pow(sqrt(m), 7.0)*dz + 1.0;
        float r = length(w);
        float b = 8 * acos(w.y / r);
        float a = 8 * atan2(w.x, w.z);
        w = p + pow(r, 8) * float3(sin(b) * sin(a), cos(b), sin(b) * cos(a));

        m = dot(w, w);
		if(m > 256.0)
            break;
    }
    return 0.25*log(m)*sqrt(m)/dz;
}

float sdDinamMandelbulb(float3 pos, float power)
{
    float3 z = pos;
    float r = 0;
    float dr = 1;
    for(int i = 0; i < 5; i++) 
    {
        r = length(z);
        if(r > 100) break;
        
        float theta = acos(z.z / r);
        float phi = atan2(z.y, z.x);
        
        dr = power * pow(r, power-1)*dr+1;
        
        r = pow(r, power);
        theta *= power;
        phi *= power;
        
        z = r * float3(sin(theta) * cos(phi), 
                sin(theta) * sin(phi), 
                cos(theta));

        z += pos;
    }
    return 0.5 * log(r) * r / dr;

}

float mandelbulb2(float3 pos, float power) {
    float3 z = pos;
	float dr = 1.0;
	float r = 0.0;
    int iterations = 0;

	for (int i = 0; i < 100 ; i++) {
        iterations = i;
		r = length(z);

		if (r>3) {
            break;
        }
        
		// convert to polar coordinates
		float theta = asin( z.z/r );
        float phi = atan2( z.y,z.x );
		dr =  pow( r, power-1.0)*power*dr + 1.0;

		// scale and rotate the point
		float zr = pow( r,power);
		theta = theta*power;
		phi = phi*power;
		
		// convert back to cartesian coordinates
		z = zr*float3( cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) );
		z+=pos;
	}
    float dst = 0.5*log(r)*r/dr;
	return dst;
}

// MengerSponge distance estimation:
// http://www.iquilezles.org/www/articles/menger/menger.htm
float MengerSponge( in float3 p )
{
   float d = sdBox(p, float3(50.0,50.0,50.0)); // size of box

   float s = 0.02; // size of squares in box
   int iterations = 0;
   
   for( int m=0; m<5; m++ )
   {
      //p = p + float3(power-1.0,power-1.0,power-1.0);
      iterations = m;
      float3 a = (p*s - 2.0 * floor(p*s/2.0))-1.0;
      s *= 3.0;
      float3 r = abs(1.0 - 3.0*abs(a));

      float da = max(r.x,r.y);
      float db = max(r.y,r.z);
      float dc = max(r.z,r.x);
      float c = (min(da,min(db,dc))-1.0)/s;

      d = max(d,c);
   }

   return d;
}

float sdJulia(float3 pos, float4 c)
{
	float4 z = float4(pos, 0);
    float md2 = 1;
    float mz2 = dot(z, z);

	[loop]
    for(int i = 0; i < 11; i++)
    {
        md2 *= 4.0 * mz2; // dz -> 2·z·dz, meaning |dz| -> 2·|z|·|dz| (can take the 4 out of the loop and do an exp2() afterwards)
        z = qsqr(z) + c; // z  -> z^2 + c

        mz2 = dot(z,z);

        if(mz2 > 4.0) break;
    }
    
    return 0.25 * sqrt(mz2/md2) * log(mz2);
}

float sdJuliabulb(float3 pos, float4 c)
{
	float3 orbit = pos;
    float dz = 1;
    
    for (int i = 0; i < 4; i++) 
    {
        float r = length(orbit);
    	float o = acos(orbit.z/r);
    	float p = atan(orbit.y/orbit.x);
        
        dz = 8*r*r*r*r*r*r*r*dz;
        
        r = r*r*r*r*r*r*r*r;
        o = 8*o;
        p = 8*p;
        
        orbit = float3(r*sin(o) * cos(p), 
                r*sin(o) * sin(p), 
                r*cos(o)) + c;
        
        if (dot(orbit, orbit) > 4.0) break;
    }
    float z = length(orbit);
    return 0.5*z*log(z)/dz;
}

float sierpinski(float3 p)
{
    const float3 va = float3(  0.0,  0.575735,  0.0 );
    const float3 vb = float3(  0.0, -1.0,  1.15470 );
    const float3 vc = float3(  1.0, -1.0, -0.57735 );
    const float3 vd = float3( -1.0, -1.0, -0.57735 );

    float a = 0;
    float s = 1;
    float r = 1;
    float dm;
    float3 v;
    [loop]
    for(int i = 0; i < 15; i++)
	{
	    float d, t;
		d = dot(p - va, p - va);

        v = va; 
        dm = d; 
        t = 0;
        
        d = dot(p - vb, p - vb); 
        if(d < dm) 
        { 
            v = vb; 
            dm=d; 
            t = 1.0; 
        }
        
        d = dot(p-vc, p-vc); 

        if(d < dm) { v = vc; dm = d; t = 2.0; }
        d = dot(p-vd,p-vd); 
        if(d < dm) { v = vd; dm = d; t = 3.0; }

		p = v + 2*(p - v); 
        r*= 2;
		a = t + 4*a; 
        s*= 4;
	}
	
	return float2((sqrt(dm)-1.0)/r, a/s);
}

float mandelbox(float3 position) {
    float SCALE = 2.75;
    float fixedRadius = 1.0;
    float FR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5;
    float MR2 = minRadius * minRadius;
    float4 scalevec = float4(SCALE, SCALE, SCALE, abs(SCALE)) / MR2;
    float C1 = abs(SCALE-1.0);
    float C2 = pow(abs(SCALE), float(1-5));
    float4 p = float4(position.xyz, 1.0); 
    float4 p0 = float4(position.xyz, 1.0);  // p.w is knighty's DEfactor
    for (int i=0; i<5; i++) {
        p.xyz = clamp(p.xyz *0.5+0.5, 0.0, 1.0) *4.0-2.0 - p.xyz; // box fold: min3, max3, mad3
        float r2 = dot(p.xyz, p.xyz);  // dp3
        p.xyzw *= clamp(max(MR2/r2, MR2), 0.0, 1.0);  // sphere fold: div1, max1.sat, mul4
        p.xyzw = p*scalevec + p0;  // mad4
    }
  return (length(p.xyz) - C1) / p.w - C2;

}

float mandelbox2(float3 p)
{
    float scale = 2;
	float3 offset = p;
	float dr = 1.0;
	for (int n = 0; n < 10; n++)
    {
		boxFold(p, dr);
		sphereFold(p, dr);
        p = scale * p + offset;
        dr = dr * abs(scale) + 1.0;
	}
	float r = length(p);
	return r / abs(dr);
}



float remap(float value, float low1, float high1, float low2, float high2)
{
    return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}

float  modc(float  a, float  b) { return a - b * floor(a/b); }
float2 modc(float2 a, float2 b) { return a - b * floor(a/b); }
float3 modc(float3 a, float3 b) { return a - b * floor(a/b); }
float4 modc(float4 a, float4 b) { return a - b * floor(a/b); }

float3 RotateX(float3 p, float angle)
{
    float rad = 0.0174532925 * angle;
    float s, c;
    sincos(rad, s, c);
    return float3(p.x, c*p.y + s*p.z, -s*p.y + c*p.z);
}
float3 RotateY(float3 v, float degree)
{
	float rad = 0.0174532925 * degree;
	float cosY = cos(rad);
	float sinY = sin(rad);
	return float3(cosY * v.x - sinY * v.z, v.y, sinY * v.x + cosY * v.z);
}
float3 RotateZ(float3 p, float angle)
{
    float rad = 0.0174532925 * angle;
    float s, c;
    sincos(rad, s, c);
    return float3(c*p.x + s*p.y, -s*p.x + c*p.y, p.z);
}

float kaleidoscopic_IFS(float3 z)
{
    int FRACT_ITER      = 20;
    float FRACT_SCALE   = 1.8;
    float FRACT_OFFSET  = 1.0;

    float c = 2.0;
    z.y = modc(z.y, c)-c/2.0;
    z = RotateZ(z, PI/2.0);
    float r;
    int n1 = 0;
    for (int n = 0; n < FRACT_ITER; n++) {
        float rotate = PI*0.5;
        z = RotateX(z, rotate);
        z = RotateY(z, rotate);
        z = RotateZ(z, rotate);

        z.xy = abs(z.xy);
        if (z.x+z.y<0.0) z.xy = -z.yx; // fold 1
        if (z.x+z.z<0.0) z.xz = -z.zx; // fold 2
        if (z.y+z.z<0.0) z.zy = -z.yz; // fold 3
        z = z*FRACT_SCALE - FRACT_OFFSET*(FRACT_SCALE-1.0);
    }
    return (length(z) ) * pow(FRACT_SCALE, -float(FRACT_ITER));
}

float2x2 rot(float a) {
	return float2x2(cos(a),sin(a),-sin(a),cos(a));	
}

float4 formula(float4 p) {
		p.xz = abs(p.xz+1.)-abs(p.xz-1.)-p.xz;
		p.y-=.25;
		float a = 35.0;
		p.x = cos(a) * p.x + sin(a) * p.y;
		p.y = -sin(a) * p.x + cos(a) * p.y;
		p=p*2./clamp(dot(p.xyz,p.xyz),.2,1.);
	return p;
}

float RemnantX(float3 pos) 
{
    float hid=0.;
	float3 tpos=pos;
	tpos.z=abs(3.-(tpos.z-6.*floor(tpos.z/6.)));
	float4 p=float4(tpos,1.);
	for (int i=0; i<4; i++) {p=formula(p);}
	float fr=(length(max(float2(0., 0.),p.yz-1.5))-1.)/p.w;
	float ro=max(abs(pos.x+1.)-.3,pos.y-.35);
		  ro=max(ro,-max(abs(pos.x+1.)-.1,pos.y-.5));
	pos.z=abs(.25-(pos.z - .5*floor(pos.z/.5)));
	
		  ro=max(ro,-max(abs(pos.z)-.2,pos.y-.3));
		  ro=max(ro,-max(abs(pos.z)-.01,-pos.y+.32));
	float d=min(fr,ro);
	return d;
}

float tglad_formula(float3 z0)
{
    z0 = modc(z0, 2.0);

    float mr=0.25, mxr=1.0;
    float4 scale=float4(-3.12,-3.12,-3.12,3.12), p0=float4(0.0,1.59,-1.0,0.0);
    float4 z = float4(z0,1.0);
    for (int n = 0; n < 3; n++) {
        z.xyz=clamp(z.xyz, -0.94, 0.94)*2.0-z.xyz;
        z*=scale/clamp(dot(z.xyz,z.xyz),mr,mxr);
        z+=p0;
    }
    float dS=(length(max(abs(z.xyz)-float3(1.2,49.0,1.4),0.0))-0.06)/z.w;
    return dS;
}


// distance function from Hartverdrahtet
// ( http://www.pouet.net/prod.php?which=59086 )
float hartverdrahtet(float3 f)
{
    float3 cs=float3(.808,.808,1.167);
    float fs=1.;
    float3 fc=0;
    float fu=10.;
    float fd=.763;
    
    // scene selection
    {
        float time = _Time.y;
        int i = int(modc(time/2.0, 9.0));
        if(i==0) cs.y=.58;
        if(i==1) cs.xy=.5;
        if(i==2) cs.xy=.5;
        if(i==3) fu=1.01,cs.x=.9;
        if(i==4) fu=1.01,cs.x=.9;
        if(i==6) cs=float3(.5,.5,1.04);
        if(i==5) fu=.9;
        if(i==7) fd=.7,fs=1.34,cs.xy=.5;
        if(i==8) fc.z=-.38;
    }
    
    //cs += sin(time)*0.2;

    float v=1.;
    for(int i=0; i<12; i++){
        f=2.*clamp(f,-cs,cs)-f;
        float c=max(fs/dot(f,f),1.);
        f*=c;
        v*=c;
        f+=fc;
    }
    float z=length(f.xy)-fu;
    return fd*max(z,abs(length(f.xy)*f.z)/sqrt(dot(f,f)))/abs(v);
}

float pseudo_kleinian(float3 p)
{
    float3 CSize = float3(0.92436,0.90756,0.92436);
    float Size = 1.0;
    float3 C = float3(0.0,0.0,0.0);
    float DEfactor=1.;
    float3 Offset = float3(0.0,0.0,0.0);
    float3 ap=p+1.;
    for(int i=0;i<10 ;i++){
        ap=p;
        p=2.*clamp(p, -CSize, CSize)-p;
        float r2 = dot(p,p);
        float k = max(Size/r2,1.);
        p *= k;
        DEfactor *= k + 0.05;
        p += C;
    }
    float r = abs(0.5*abs(p.z-Offset.z)/DEfactor);
    return r;
}

float pseudo_knightyan(float3 p)
{
    float3 CSize = float3(0.63248,0.78632,0.875);
    float DEfactor=1.;
    for(int i=0;i<6;i++){
        p = 2.*clamp(p, -CSize, CSize)-p;
        float k = max(0.70968/dot(p,p),1.);
        p *= k;
        DEfactor *= k + 0.05;
    }
    float rxy=length(p.xy);
    return max(rxy-0.92784, abs(rxy*p.z) / length(p))/DEfactor;
}

sampler2D g_depth_prev;
sampler2D g_depth;
sampler2D g_velocity;

float cross_depth_sample(float2 t, sampler2D s, float o)
{
    float2 p = (_ScreenParams.zw - 1.0)*o;
    float d1 = tex2D(s, t).x;
    float d2 = min(
        min(tex2D(s, t+float2( p.x, 0.0)).x, tex2D(s, t+float2(-p.x, 0.0))).x,
        min(tex2D(s, t+float2( 0.0, p.y)).x, tex2D(s, t+float2( 0.0,-p.y))).x );
    return min(d1, d2);
}

float sample_prev_depth(float2 t)
{
    return max(tex2D(g_depth_prev, t).x-0.001, _ProjectionParams.y);
}

float sample_upper_depth(float2 t)
{
    return max(cross_depth_sample(t, g_depth, 2.0)*0.995, _ProjectionParams.y);
}

float  iq_rand(float  p)
{
    return frac(sin(p)*43758.5453);
}
float2 iq_rand(float2 p)
{
    p = float2(dot(p, float2(127.1, 311.7)), dot(p, float2(269.5, 183.3)));
    return frac(sin(p)*43758.5453);
}
float3 iq_rand(float3 p)
{
    p = float3(dot(p, float3(127.1, 311.7, 311.7)), dot(p, float3(269.5, 183.3, 183.3)), dot(p, float3(269.5, 183.3, 183.3)));
    return frac(sin(p)*43758.5453);
}


float hash(float2 p)
{
    float h = dot(p, float2(127.1, 311.7));
    return frac(sin(h)*43758.5453123);
}
float hash(float3 p)
{
    float h = dot(p, float3(127.1, 311.7, 496.3));
    return frac(sin(h)*43758.5453123);
}

float iqnoise(float2 p)
{
    float2 i = floor(p);
    float2 f = frac(p);
    float2 u = f*f*(3.0 - 2.0*f);
    return -1.0 + 2.0*lerp(lerp(hash(i + float2(0.0, 0.0)),
        hash(i + float2(1.0, 0.0)), u.x),
        lerp(hash(i + float2(0.0, 1.0)),
        hash(i + float2(1.0, 1.0)), u.x), u.y);
}

float iqnoise(float3 p)
{
    float3 i = floor(p);
    float3 f = frac(p);
    float3 u = f*f*(3.0 - 2.0*f);
    return -1.0 + 2.0*lerp(lerp(hash(i + float2(0.0, 0.0)),
        hash(i + float2(1.0, 0.0)), u.x),
        lerp(hash(i + float2(0.0, 1.0)),
        hash(i + float2(1.0, 1.0)), u.x), u.y);
}


// trinoise

float tri(float x)
{
    return abs(frac(x) - .5);
}

float3 tri3(float3 p)
{
    return float3(
        tri(p.z + tri(p.y * 1.)),
        tri(p.z + tri(p.x * 1.)),
        tri(p.y + tri(p.x * 1.)) );
}

float trinoise(float3 p, float spd, float time)
{
    float z = 1.4;
    float rz = 0.;
    float3  bp = p;
    for (float i = 0.; i <= 3.; i++) {
        float3 dg = tri3(bp * 2.);
        p += (dg + time * .1 * spd);
        bp *= 1.8;
        z *= 1.5;
        p *= 1.2;
        float t = tri(p.z + tri(p.x + tri(p.y)));
        rz += t / z;
        bp += 0.14;
    }
    return rz;
}

float trinoise(float3 p)
{
    return trinoise(p, 0.0, 0.0);
}
