using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using Vector3 = UnityEngine.Vector3;

[ExecuteInEditMode]
public class RaymarchCollider : MonoBehaviour
{
    public float de;
    public Vector3 position;
    SphereCollider myCollider;
    Rigidbody rb;
    
    private void Start()
    {
        myCollider = GetComponent<SphereCollider>();
        rb = GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    void Update()
    {
        position = GetComponent<Transform>().position;
        // Okay so, there is progress here, but the distance function thinks the sponge fractal is a solid cube, and doesnt take into account the intricacies of the fractal (Maybe its with the sphere object itself??)
        // Should try and test it out on other fractals and see what happens
        de = this.MengerSponge(position);;
        //Check if the distance estimate indicates a collision
        if (de <= 25)
        {
            // Push object back to surface of fractal 
            rb.AddForce(transform.up * 50f);
        }
    }
    
    float sdBox(Vector3 p, Vector3 b)
    {
        p.x = Math.Abs(p.x);
        p.y = Math.Abs(p.y);
        p.z = Math.Abs(p.z);
        Vector3 d = p - b;

        Vector3 e;
        e.x = Math.Max(d.x, 0);
        e.y = Math.Max(d.y, 0);
        e.z = Math.Max(d.z, 0);

        float length = (float) Math.Sqrt((e.x * e.x) + (e.y * e.y) + (e.z * e.z));
        
        return (float)(Math.Min(Math.Max(d.x, Math.Max(d.y, d.z)), 0.0) +
              length);
    }
    
    float MengerSponge( Vector3 p )
    {
        float d = sdBox(p, new Vector3(50,50,50)); // size of box

        float s = 0.02f; // size of squares in box

        for( int m=0; m<5; m++ )
        {
            Vector3 a;
            a.x = (float)((p.x * s - 2.0 * Math.Floor(p.x * s / 2.0f)) - 1.0);
            a.y = (float)((p.y * s - 2.0 * Math.Floor(p.y * s / 2.0f)) - 1.0);
            a.z = (float)((p.z * s - 2.0 * Math.Floor(p.z * s / 2.0f)) - 1.0);
            
            s *= 3.0f;

            Vector3 r;

            r.x = Math.Abs(1 - 3f * Math.Abs(a.x));
            r.y = Math.Abs(1 - 3f * Math.Abs(a.y));
            r.z = Math.Abs(1 - 3f * Math.Abs(a.z));

            float da = Math.Max(r.x,r.y);
            float db = Math.Max(r.y,r.z);
            float dc = Math.Max(r.z,r.x);
            float c = (float)(Math.Min(da,Math.Max(db,dc))-1.0)/s;

            d = Math.Max(d,c);
        }

        return d;
    }
}
