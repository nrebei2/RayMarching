Shader "Noah/RaymarchShader"
{
	Properties
	{
		_MainTex ("Texture", 2D) = "white" {}
	}
	SubShader
	{
		// No culling or depth
		Cull Off ZWrite Off ZTest Always

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma target 3.0
			#define PI 3.14159265358979323846
			
			#include "UnityCG.cginc"
			#include "DistanceFunctions.cginc"

			sampler2D _MainTex;
			uniform sampler2D _CameraDepthTexture;
			uniform float4x4 _CamFrustum, _CamToWorld;
			uniform int _MaxIterations;
			uniform float _Accuracy;
			uniform float _maxDistance, _box1round, _boxSphereSmooth, _sphereIntersectSmooth;
			uniform float4 _sphere1, _sphere2, _box1;
			
			uniform float3 _LightDir, _LightCol;
			uniform float _LightIntensity;
			
			uniform fixed4 _FractalColor;
			
			uniform fixed4 _GroundColor;
			uniform float _ColorIntensity;
			
			uniform float2 _ShadowDistance;
			uniform float _ShadowIntensity, _ShadowPenumbra;

			uniform int _ReflectionCount;
			uniform float _ReflectionIntensity;
			uniform float _EnvReflIntenisty;
			uniform samplerCUBE _ReflectionCube;
			
			uniform float4 _fractal;
			uniform float _fractalSmooth;
			uniform float3 _fractaldegreeRotate;
			
			uniform float3 _modInterval;
			
			uniform float3 _modBool;
			
			uniform float _fractalNumber;
            
            uniform float _power;
            

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float2 uv : TEXCOORD0;
				float4 vertex : SV_POSITION;
				float3 ray : TEXCOORD1;
			};

			v2f vert (appdata v)
			{
				v2f o;
				half index = v.vertex.z;
				v.vertex.z = 0;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = v.uv;

				o.ray = _CamFrustum[(int)index].xyz;

				o.ray /= abs(o.ray.z);

				o.ray = mul(_CamToWorld, o.ray);

				return o;
			}
            
			//float BoxSphere(float3 p)
			//{
			//	float Sphere1 = sdSphere(p - _sphere1.xyz, _sphere1.w);
			//	float Box1 = sdRoundBox(p - _box1.xyz, _box1.www, _box1round);
			//	float combine1 = opSS(Sphere1, Box1, _boxSphereSmooth);
			//	float Sphere2 = sdSphere(p - _sphere2.xyz, _sphere2.w);
			//	float combine2 = opIS(Sphere2, combine1, _sphereIntersectSmooth);
			//
			//	return combine2;
			//}

			

			float4 distanceField(float3 p)
			{
			    p = RotateZ(RotateY(RotateX(p - _fractal.xyz, _fractaldegreeRotate.x), _fractaldegreeRotate.y), _fractaldegreeRotate.z);

			    if (_modBool.x == 1) {
    			    float modX = pMod1(p.x, _modInterval.x);
			    }
			    if (_modBool.y == 1) {
			        float modY = pMod1(p.y, _modInterval.y);
			    }
			    if (_modBool.z == 1) {
			        float modZ = pMod1(p.z, _modInterval.z);
			    }
			    
			    float4 fractal;
			    float4 c = 0.45* cos( float4(0.5,3.9,1.4,1.1) + _Time.y * float4(1.2,1.7,1.3,2.5) ) - float4(0.3,0.0,0.0,0.0);
			    float sinTime = sin(_Time.y / 1);
				float power = remap(sinTime, -1, 1, 4, 9);
			    
			    switch (_fractalNumber)
                {
                    case 0:
                        fractal = float4(_FractalColor.rgb, sdMandelbulb(p));
                        break;
                    case 1:
                        fractal = float4(_FractalColor.rgb, sdDinamMandelbulb(p, power));
                        break;
                    case 2:
                        fractal = float4(_FractalColor.rgb, sdJulia(p, c));
                        break;
                    case 3:
                        fractal = float4(_FractalColor.rgb, sdJuliabulb(p, c));
                        break;
                    case 4:
                        fractal = float4(_FractalColor.rgb, sierpinski(p));
                        break;
                    case 5:
                        fractal = float4(_FractalColor.rgb, mandelbox(p));
                        break;
                    case 6:
                        fractal = float4(_FractalColor.rgb, kaleidoscopic_IFS(p));
                        break;
                    case 7:
                        fractal = float4(_FractalColor.rgb, tglad_formula(p));
                        break;
                    case 8:
                        fractal = float4(_FractalColor.rgb, hartverdrahtet(p));
                        break;
                    case 9:
                        fractal = float4(_FractalColor.rgb, pseudo_kleinian(p));
                        break;
                    case 10:
                        fractal = float4(_FractalColor.rgb, pseudo_knightyan(p));
                        break;
                    case 11:
                        fractal = float4(_FractalColor.rgb, mandelbulb2(p, _power));
                        break;        
                    case 12:
                        fractal = float4(_FractalColor.rgb, MengerSponge(p));
                        break;      
                    case 13:
                        fractal = float4(_FractalColor.rgb, apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.)));
                        break; 
                    case 14:
                        fractal = float4(_FractalColor.rgb, sdPlane(p, float4(0,1,0,0)));
                        break;    
                    case 15:
                        fractal = float4(_FractalColor.rgb, FCT_BBSK(p));     
                        break;                          
                }
			    
                return fractal;
			}

			float3 getNormal(float3 p)
			{
				const float2 offset = float2(0.001, 0.0);
				float3 n = float3(
					distanceField(p + offset.xyy).w - distanceField(p - offset.xyy).w,
					distanceField(p + offset.yxy).w - distanceField(p - offset.yxy).w,
					distanceField(p + offset.yyx).w - distanceField(p - offset.yyx).w);
				return normalize(n);
			}

			float hardShadow(float3 ro, float rd, float mint, float maxt)
			{
				for (float t = mint; t < maxt;)
				{
					float h = distanceField(ro + rd * t).w;
					if (h < 0.001)
					{
						return 0.0;
					}
					t += h;
				}
				return 1.0;
			}
			float softShadow(float3 ro, float rd, float mint, float maxt, float k)
			{
				float result = 1.0;
				for (float t = mint; t < maxt;)
				{
					float h = distanceField(ro + rd * t).w;
					if (h < 0.001)
					{
						return 0.0;
					}
					result = min(result, k*h / t);
					t += h;
				}
				return result;
			}

			uniform float _AoStepsize, _AoInstensity;
			uniform int _AoIterations;

			float AmbientOcclusion(float3 p, float3 n)
			{
				float step = _AoStepsize;
				float ao = 0.0;
				float dist;
				for (int i = 1; i <= _AoIterations; i++)
				{
					dist = step * i;
					ao += max(0.0, (dist - distanceField(p + n * dist).w) / dist);
				}
				return (1.0 - ao * _AoInstensity);
			}
			
			float3 Shading(float3 p, float3 n, fixed3 c)
			{
				float3 result;
				//Diffuse Color
				float3 color = c.rgb * _ColorIntensity;
				//Directional Light
				float3 light = (_LightCol * dot(-_LightDir, n) * 0.5 + 0.5) * _LightIntensity;
				//Shadows
				float shadow = softShadow(p, -_LightDir, _ShadowDistance.x, _ShadowDistance.y, _ShadowPenumbra) * 0.5 + 0.5;
				shadow = max(0.0, pow(shadow, _ShadowIntensity));
				//AMBIENNT OOCCLLUUU
				float ao = AmbientOcclusion(p, n);


				result = color * light * shadow * ao;

				return result;
			}
			
			bool raymarching(float3 ro, float3 rd, float depth, float maxDistance, int maxIterations, inout float3 p, inout fixed3 dColor)
			{
				bool hit;

				float t = 0; //distance travelled along the ray direction

				for (int i = 0; i < maxIterations; i++)
				{
					if (t > maxDistance || t >= depth)
					{
						//env
						hit = false;
						break;
					}

					p = ro + rd * t;
					//check for hit in distfield
					float4 d = distanceField(p);

					if (d.w < _Accuracy)
					{
					    dColor = d.rgb;
						hit = true;
						break;
					}
					t += d.w;
				}

				return hit;
			}


			fixed4 frag (v2f i) : SV_Target
			{
				float depth = LinearEyeDepth(tex2D(_CameraDepthTexture, i.uv).r);
				depth *= length(i.ray);
				fixed3 col = tex2D(_MainTex, i.uv);
				float3 rayDirection = normalize(i.ray.xyz);
				float3 rayOrigin = _WorldSpaceCameraPos;
				fixed4 result;
				float3 hitPosition;
				fixed3 dColor;

				bool hit = raymarching(rayOrigin, rayDirection, depth, _maxDistance, _MaxIterations, hitPosition, dColor);
				if (hit)
				{				
					//shading!
					float3 n = getNormal(hitPosition);
					float3 s = Shading(hitPosition, n, dColor);

					result = fixed4(s, 1);
					result += fixed4(texCUBE(_ReflectionCube, n).rgb * _EnvReflIntenisty *_ReflectionIntensity, 0);
					
					if (_ReflectionCount > 0)
					{
						rayDirection = normalize(reflect(rayDirection, n));
						rayOrigin = hitPosition + (rayDirection * 0.01);
                        hit = raymarching(rayOrigin, rayDirection, _maxDistance, _maxDistance * 0.5, _MaxIterations / 2.0, hitPosition, dColor);
                        if (hit) {
                            float3 n = getNormal(hitPosition);
					        float3 s = Shading(hitPosition, n, dColor);
					        result += fixed4(s * _ReflectionIntensity, 0);
					        if (_ReflectionCount > 1) {
					            rayDirection = normalize(reflect(rayDirection, n));
						        rayOrigin = hitPosition + (rayDirection * 0.01);
						        hit = raymarching(rayOrigin, rayDirection, _maxDistance, _maxDistance * 0.25, _MaxIterations / 4.0, hitPosition, dColor);
						        if (hit) {
						            float3 n = getNormal(hitPosition);
					                float3 s = Shading(hitPosition, n, dColor);
					                result += fixed4(s * _ReflectionIntensity * 0.5, 0);
						        }
					        }
                        }
					}
				}
				else //miss
				{
					result = fixed4(0, 0, 0, 0);
				}


				return fixed4(col * (1.0 - result.w) + result.xyz * result.w,1.0);
			}
			ENDCG
		}
	}
}
