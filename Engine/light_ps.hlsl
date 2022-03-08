// Light pixel shader
// Calculate diffuse lighting for a single directional light(also texturing)

Texture2D shaderTexture : register(t0);
Texture2D shaderTexture2 : register(t1);
SamplerState SampleType : register(s0);

static const float SNOW = 0.0000002f;


cbuffer LightBuffer : register(b0)
{
	float4 ambientColor;
    float4 diffuseColor;
    float3 lightPosition;
    float padding;
};

struct InputType
{
    float4 position : SV_POSITION;
    float2 tex : TEXCOORD0;
    float3 normal : NORMAL;
	float3 position3D : TEXCOORD2;
	float height : FLOAT;
};

float4 main(InputType input) : SV_TARGET
{
	float4	textureColor;
float4	textureColor2;
    float3	lightDir;
    float	lightIntensity;
    float4	color;

	// Invert the light direction for calculations.
	lightDir = normalize(input.position3D - lightPosition);

	// Calculate the amount of light on this pixel.
	lightIntensity = saturate(dot(input.normal, -lightDir));

	// Determine the final amount of diffuse color based on the diffuse color combined with the light intensity.
	color = ambientColor + (diffuseColor * lightIntensity); //adding ambient
	color = saturate(color);

	if (input.height > SNOW)
	{
		textureColor = shaderTexture2.Sample(SampleType, input.tex);
	}
	else if(input.height > 0.000005f && input.height < 0.000002f){

		textureColor = shaderTexture.Sample(SampleType, input.tex) / 2;
		textureColor2 = shaderTexture2.Sample(SampleType, input.tex) / 2;
		textureColor = textureColor + textureColor2;

	}
	else {
		// Sample the pixel color from the texture using the sampler at this texture coordinate location.
		textureColor = shaderTexture.Sample(SampleType, input.tex);
	}

	color = color * textureColor;

    return color;
}

