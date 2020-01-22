using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(Camera))]
[ExecuteInEditMode]
public class RaymarchingCamera : SceneViewFilter
{

    [SerializeField] private Shader _shader;

    public Material _raymarchMaterial
    {
        get
        {
            if (!_raymarchMat && _shader)
            {
                _raymarchMat = new Material(_shader);
                _raymarchMat.hideFlags = HideFlags.HideAndDontSave;
            }
            return _raymarchMat;
        }
    }

    private Material _raymarchMat;
    
    public Camera _camera
    {
        get
        {
            if(!_cam)
            {
                _cam = GetComponent<Camera>();
            }
            return _cam;
        }
    }

    private Camera _cam;
    [Header("Setup")]
    public float _maxDistance;
    [Range(1,3000)]
    public int _MaxIterations;
    [Range(0.1f, 0.0001f)]
    public float _Accuracy;

    [Header("Directional Light")]
    public Transform _directionalLight;
    public Color _LightCol;
    public float _LightIntensity;


    [Header("Shading")]
    [Range(0, 4)]
    public float _ShadowIntensity;
    public Vector2 _ShadowDistance;
    [Range(1,128)]
    public float _ShadowPenumbra;

    [Header("Ambient Occlusion")]
    [Range(0.01f, 10.0f)]
    public float _AoStepsize;
    [Range(1,5)]
    public int _AoIterations;
    [Range(0,1)]
    public float _AoInstensity;

    [Header("Reflection")]
    [Range(0, 2)]
    public int _ReflectionCount;
    [Range(0, 1)]
    public float _ReflectionIntensity;
    [Range(0, 1)]
    public float _EnvReflIntenisty;

    public Cubemap _ReflectionCube;

    public enum Fractal{ StaticMandelbulb, DynamicMandelbulb, Julia, Juliabulb, Sierpinski, Mandelbox, 
        KaleidoscopicIFS, Tglad, Hartverdrahtet, PseudoKleinian, PseudoKnightyan, Mandelbulb2, MengerSponge, apo }
    
    [Header("Signed Distance Field")]
    public Fractal fractal;
    private float _fractalNumber;
    public Vector4 _fractal;
    public float _fractalSmooth;
    public Vector3 _fractaldegreeRotate;
    public float _power;
    

    [Header("Color")] 
    public Color _FractalColor;
    [Range(0, 4)] 
    public float _ColorIntensity;
    
    [Header("Make Object Repeat Indefinitely on Axis (Value 1 is True, Any Other Value is False)")]
    public Vector3 _modBool;
    
    [Header("ModInterval")]
    public Vector3 _modInterval;

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        if(!_raymarchMaterial)
        {
            Graphics.Blit(source, destination);
            return;
        }

        switch (fractal)
        {
            case Fractal.StaticMandelbulb:
                _fractalNumber = 0;
                break;
            case Fractal.DynamicMandelbulb:
                _fractalNumber = 1;
                break;
            case Fractal.Julia:
                _fractalNumber = 2;
                break;
            case Fractal.Juliabulb:
                _fractalNumber = 3;
                break;
            case Fractal.Sierpinski:
                _fractalNumber = 4;
                break;
            case Fractal.Mandelbox:
                _fractalNumber = 5;
                break;
            case Fractal.KaleidoscopicIFS:
                _fractalNumber = 6;
                break;
            case Fractal.Tglad:
                _fractalNumber = 7;
                break;
            case Fractal.Hartverdrahtet:
                _fractalNumber = 8;
                break;
            case Fractal.PseudoKleinian:
                _fractalNumber = 9;
                break;
            case Fractal.PseudoKnightyan:
                _fractalNumber = 10;
                break;
            case Fractal.Mandelbulb2:
                _fractalNumber = 11;
                break;
            case Fractal.MengerSponge:
                _fractalNumber = 12;
                break;
            case Fractal.apo:
                _fractalNumber = 13;
                break;
        }

        _raymarchMaterial.SetFloat("_fractalNumber", _fractalNumber);

        _raymarchMaterial.SetVector("_modBool", _modBool);
        _raymarchMaterial.SetVector("_modInterval", _modInterval);

        _raymarchMaterial.SetColor("_FractalColor", _FractalColor);
        _raymarchMaterial.SetFloat("_ColorIntensity", _ColorIntensity);
        
        _raymarchMaterial.SetVector("_LightDir", _directionalLight ? _directionalLight.forward : Vector3.down);
        _raymarchMaterial.SetColor("_LightCol", _LightCol);
        _raymarchMaterial.SetFloat("_LightIntensity", _LightIntensity);
        _raymarchMaterial.SetFloat("_ShadowIntensity", _ShadowIntensity);
        _raymarchMaterial.SetFloat("_ShadowPenumbra", _ShadowPenumbra);
        _raymarchMaterial.SetVector("_ShadowDistance", _ShadowDistance);
        _raymarchMaterial.SetMatrix("_CamFrustum", CamFrustrum(_camera));
        _raymarchMaterial.SetMatrix("_CamToWorld", _camera.cameraToWorldMatrix);
        _raymarchMaterial.SetFloat("_maxDistance", _maxDistance);
        _raymarchMaterial.SetFloat("_Accuracy", _Accuracy);
        _raymarchMaterial.SetInt("_MaxIterations", _MaxIterations);

        _raymarchMaterial.SetVector("_fractal", _fractal);
        _raymarchMaterial.SetFloat("_fractalSmooth", _fractalSmooth);
        _raymarchMaterial.SetVector("_fractaldegreeRotate", _fractaldegreeRotate);
        _raymarchMaterial.SetFloat("_power", _power);

        _raymarchMaterial.SetFloat("_AoStepsize", _AoStepsize);
        _raymarchMaterial.SetFloat("_AoInstensity", _AoInstensity);
        _raymarchMaterial.SetInt("_AoIterations", _AoIterations);

        _raymarchMaterial.SetInt("_ReflectionCount", _ReflectionCount);
        _raymarchMaterial.SetFloat("_ReflectionIntensity", _ReflectionIntensity);
        _raymarchMaterial.SetFloat("_EnvReflIntenisty", _EnvReflIntenisty);
        _raymarchMaterial.SetTexture("_ReflectionCube", _ReflectionCube);

        RenderTexture.active = destination;
        _raymarchMaterial.SetTexture("_MainTex", source);
        GL.PushMatrix();
        GL.LoadOrtho();
        _raymarchMaterial.SetPass(0);
        GL.Begin(GL.QUADS);

        //BL
        GL.MultiTexCoord2(0, 0.0f, 0.0f);
        GL.Vertex3(0.0f, 0.0f, 3.0f);
        //BR
        GL.MultiTexCoord2(0, 1.0f, 0.0f);
        GL.Vertex3(1.0f, 0.0f, 2.0f);
        //TR
        GL.MultiTexCoord2(0, 1.0f, 1.0f);
        GL.Vertex3(1.0f, 1.0f, 1.0f);
        //TL
        GL.MultiTexCoord2(0, 0.0f, 1.0f);
        GL.Vertex3(0.0f, 1.0f, 0.0f);

        GL.End();
        GL.PopMatrix();
    }

    private Matrix4x4 CamFrustrum(Camera cam)
    {
        Matrix4x4 frustum = Matrix4x4.identity;
        float fov = Mathf.Tan((cam.fieldOfView * 0.5f) * Mathf.Deg2Rad);

        Vector3 goUp = Vector3.up * fov;
        Vector3 goRight = Vector3.right * fov * cam.aspect;

        Vector3 TL = (-Vector3.forward - goRight + goUp);
        Vector3 TR = (-Vector3.forward + goRight + goUp);
        Vector3 BR = (-Vector3.forward + goRight - goUp);
        Vector3 BL = (-Vector3.forward - goRight - goUp);

        frustum.SetRow(0, TL);
        frustum.SetRow(1, TR);
        frustum.SetRow(2, BR);
        frustum.SetRow(3, BL);

        return frustum;
    }
}
