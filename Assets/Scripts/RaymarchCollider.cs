using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;
using Vector4 = UnityEngine.Vector4;

[ExecuteInEditMode]
public class RaymarchCollider : MonoBehaviour
{
    public float de;
    public Vector3 position;
    SphereCollider myCollider;
    Rigidbody rb;
    public Vector3 sdfNormal;
    public bool okToJump;
    
    Vector3 _gravityVec = new Vector3(0,-0.7071f,0.7071f);
    
    // TODO: do definitions of shader code and turn it to readable C# code, making it easy to translate code from shadercode to c#
    // TODO: implement collision similar to Yeomada
    // TODO: implement the sdf FCT_BBSK from distancefuction to see if something is really wrong with what im doing (try this first)
    private void Start()
    {
        myCollider = GetComponent<SphereCollider>();
        rb = GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    void Update()
    {        
        position = GetComponent<Transform>().position;
        okToJump = false;
        Vector3 playerV = rb.velocity;
        //sdfNormal = getNormal(position);

        // Okay so, there is progress here, but the distance function thinks the sponge fractal is a solid cube, and doesnt take into account the intricacies of the fractal (Maybe its with the sphere object itself??)
        // Should try and test it out on other fractals and see what happens
        de = FCT_BBSK(position);;
        //Check if the distance estimate indicates a collision
        if (de <= 0.5)
        {
            // Push object back to surface of fractal
            okToJump = true;
            rb.AddForce(playerV * -1 , ForceMode.Impulse);
            //rb.AddForce(FeetN*(0.3-FeetDist)*rb.mass * 16 * invAbsFdot);
        }
    }

    float sdPlane(Vector3 p, Vector4 n) {
        // n must be normalized
        // n is the normal, meaning n = float(0,1,0,0) is pointing up on the y
        return Vector3.Dot(p, new Vector3(n.x, n.y, n.z)) + n.w;
    }
    
    // TODO: fix this mess
    /*
    Vector3 getNormal(Vector3 p)
    {
        
        Vector2 offset = new Vector2(0.001f, 0.0f);
        Vector3 n = new Vector3(
            FCT_BBSK(new Vector3(p.x + offset.x, p.y + offset.y, p.z + offset.y)) - FCT_BBSK(p - offset.xyy).w,
            FCT_BBSK(p + offset.yxy).w - FCT_BBSK(p - offset.yxy).w,
            FCT_BBSK(p + offset.yyx).w - FCT_BBSK(p - offset.yyx).w);
        return Vector3.Normalize(n);
    }
    */
    
    float sdBox(Vector3 p, Vector3 b)
    {
        p.x = Mathf.Abs(p.x);
        p.y = Mathf.Abs(p.y);
        p.z = Mathf.Abs(p.z);
        Vector3 d = p - b;

        Vector3 e;
        e.x = Mathf.Max(d.x, 0);
        e.y = Mathf.Max(d.y, 0);
        e.z = Mathf.Max(d.z, 0);

        float length = (float) Mathf.Sqrt((e.x * e.x) + (e.y * e.y) + (e.z * e.z));
        
        return (float)(Mathf.Min(Mathf.Max(d.x, Mathf.Max(d.y, d.z)), 0) +
                       length);
    }
    
    float MengerSponge( Vector3 p )
    {
        float d = sdBox(p, new Vector3(50,50,50)); // size of box

        float s = 0.02f; // size of squares in box

        for( int m=0; m<5; m++ )
        {
            Vector3 a;
            a.x = (float)((p.x * s - 2.0 * Mathf.Floor(p.x * s / 2.0f)) - 1.0);
            a.y = (float)((p.y * s - 2.0 * Mathf.Floor(p.y * s / 2.0f)) - 1.0);
            a.z = (float)((p.z * s - 2.0 * Mathf.Floor(p.z * s / 2.0f)) - 1.0);
            
            s *= 3.0f;

            Vector3 r;

            r.x = Mathf.Abs(1 - 3f * Mathf.Abs(a.x));
            r.y = Mathf.Abs(1 - 3f * Mathf.Abs(a.y));
            r.z = Mathf.Abs(1 - 3f * Mathf.Abs(a.z));

            float da = Mathf.Max(r.x,r.y);
            float db = Mathf.Max(r.y,r.z);
            float dc = Mathf.Max(r.z,r.x);
            float c = (float)(Mathf.Min(da,Mathf.Max(db,dc))-1.0)/s;

            d = Mathf.Max(d,c);
        }

        return d;
    }
    
    float FCT_BBSK(Vector3 pos) {
        Vector3 cFcParams = new Vector3(2.18f, -0.18f, 0);
        Vector3 CSize = new Vector3(1.4f,0.87f, 1.1f);
        Vector3 p;
        p.x = 2 * pos.x;
        p.y = 2 * pos.z;
        p.z = 2 * pos.y;
        float scale = 1.0f;
    
        for( int i=0; i < 4;i++ ) 
        {
            p.x = 2 * Mathf.Clamp(p.x, -CSize.x, CSize.x) - p.x;
            p.y = 2 * Mathf.Clamp(p.y, -CSize.y, CSize.y) - p.y;
            p.z = 2 * Mathf.Clamp(p.z, -CSize.z, CSize.x) - p.z;
            float r2 = Vector3.Dot(p,p);
            //float r2 = dot(p,p+sin(p.z*.5)); //Alternate fractal
            float k = Mathf.Max((2f)/(r2), .17f);
            p *= k;
            //p *=rot;
            //p= p.yzx;
            p+=new Vector3(0.2f,0.2f,-0.5f);
            scale *= k;
        }
        
        p.x = 2 * Mathf.Clamp(p.x, -CSize.x * 4, CSize.x * 4) - p.x;
        p.y = 2 * Mathf.Clamp(p.y, -CSize.y * 4, CSize.y * 4) - p.y;
        p.z = 2 * Mathf.Clamp(p.z, -CSize.z * 4, CSize.x * 4) - p.z;

        for( int i=0; i < 8;i++ )
        {
            p.x = 2 * Mathf.Clamp(p.x, -CSize.x, CSize.x) - p.x;
            p.y = 2 * Mathf.Clamp(p.y, -CSize.y, CSize.y) - p.y;
            p.z = 2 * Mathf.Clamp(p.z, -CSize.z, CSize.x) - p.z;
            float r2 = Vector3.Dot(p,p);
            //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
            float k = Math.Max((cFcParams.x)/(r2),  0.027f);
            p     *= k;
            scale *= k;
            p.y += cFcParams.y;
        }

        float l = Mathf.Sqrt(Vector2.SqrMagnitude(new Vector2(p.x, p.y)));
        //l = mix(l,l2,0.5);
        float rxy = l - 4;
        float n = p.z;
        rxy = Math.Max(rxy, -(n) / 4);
        float dist = (rxy) / Math.Abs(scale);
        dist *= .75f;

        return dist;

    }


}
