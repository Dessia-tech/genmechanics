from string import Template

babylon_template = Template(
"""
<!doctype html>
<html>
<head>
   <meta charset="utf-8">
   <title>$name</title>
   <style>
      html, body {
         overflow: hidden;
         width: 100%;
         height: 100%;
         margin: 0;
         padding: 0;
      }
      #renderCanvas {
         width: 100%;
         height: 100%;
         touch-action: none;
      }
   </style>
    <!-- Link to the last version of BabylonJS -->

      <script src="https://preview.babylonjs.com/babylon.js"></script>
      <script src="https://preview.babylonjs.com/loaders/babylonjs.loaders.min.js"></script>
      <script src="https://code.jquery.com/pep/0.4.3/pep.js"></script>
      <script src='https://unpkg.com/earcut@2.1.1/dist/earcut.min.js'></script>
      <script src='https://preview.babylonjs.com/gui/babylon.gui.min.js'></script>


</head>
<body>
   <canvas id="renderCanvas"></canvas>
   <script type="text/javascript">
      // Get the canvas element from our HTML below
      var canvas = document.querySelector("#renderCanvas");
      // Load the BABYLON 3D engine
      var engine = new BABYLON.Engine(canvas, true);
      // -------------------------------------------------------------
      // Here begins a function that we will 'call' just after it's built

var createScene = function () {
    // This creates a basic Babylon Scene object (non-mesh)
        var scene = new BABYLON.Scene(engine);
        scene.useRightHandedSystem = true;
	var cam = new BABYLON.ArcRotateCamera("ArcRotateCamera", 0, 0, 2*$length, new BABYLON.Vector3($center), scene);
	cam.wheelPrecision=50./$length;
	cam.pinchPrecision=50./$length;
	cam.panningSensibility=250./$length;
	cam.minZ=0.01*$length;
	cam.attachControl(canvas);
	cam.inertia = 0;
	cam.panningInertia = 0;
  cam.upVector = new BABYLON.Vector3(0, 0, 1);
	// cam.mode = BABYLON.Camera.ORTHOGRAPHIC_CAMERA;

  var red_material = new BABYLON.StandardMaterial("red_material", scene);
  red_material.diffuseColor = new BABYLON.Color3(0.8, 0, 0);
  red_material.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);
  red_material.emissiveColor = new BABYLON.Color3(0.8, 0, 0);
  red_material.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);

  var green_material = new BABYLON.StandardMaterial("red_material", scene);
  green_material.diffuseColor = new BABYLON.Color3(0, 0.8, 0);
  green_material.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);
  green_material.emissiveColor = new BABYLON.Color3(0, 0.8, 0);
  green_material.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);

  var blue_material = new BABYLON.StandardMaterial("red_material", scene);
  blue_material.diffuseColor = new BABYLON.Color3(0, 0, 0.8);
  blue_material.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);
  blue_material.emissiveColor = new BABYLON.Color3(0, 0, 0.8);
  blue_material.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);


	var light1 = new BABYLON.HemisphericLight("light1", new BABYLON.Vector3(-1, -1, -1), scene);
	light1.intensity=0.5;
	light1.specular = new BABYLON.Color3(0, 0, 0);

  var light2 = new BABYLON.PointLight("light2", new BABYLON.Vector3(0, 0, 0), scene);
  light2.specular = new BABYLON.Color3(0, 0, 0);
  light2.intensity = 0.3;
  light2.parent = cam;

  var light3 = new BABYLON.HemisphericLight("light3", new BABYLON.Vector3(1, 1, 1), scene);
  light3.specular = new BABYLON.Color3(0, 0, 0);
  light3.intensity = 0.50;


    var showAxis = function (size) {
        var makeTextPlane = function (text, color, size) {
            var dynamicTexture = new BABYLON.DynamicTexture("DynamicTexture", 50, scene, true);
            dynamicTexture.hasAlpha = true;
            dynamicTexture.drawText(text, 5, 40, "bold 36px Arial", color, "transparent", true);
            var plane = new BABYLON.Mesh.CreatePlane("TextPlane", size, scene, true);
            plane.material = new BABYLON.StandardMaterial("TextPlaneMaterial", scene);
            plane.material.backFaceCulling = false;
            plane.material.specularColor = new BABYLON.Color3(0, 0, 0);
            plane.material.diffuseTexture = dynamicTexture;
            return plane;
        };

        var axisX = BABYLON.Mesh.CreateLines("axisX", [
            new BABYLON.Vector3.Zero(), new BABYLON.Vector3(size, 0, 0), new BABYLON.Vector3(size * 0.95, 0.05 * size, 0),
            new BABYLON.Vector3(size, 0, 0), new BABYLON.Vector3(size * 0.95, -0.05 * size, 0)
        ], scene);
        axisX.color = new BABYLON.Color3(1, 0, 0);
        var xChar = makeTextPlane("X", "red", size / 10);
        xChar.position = new BABYLON.Vector3(0.9 * size, -0.05 * size, 0);

        var axisY = BABYLON.Mesh.CreateLines("axisY", [
            new BABYLON.Vector3.Zero(), new BABYLON.Vector3(0, size, 0), new BABYLON.Vector3(-0.05 * size, size * 0.95, 0),
            new BABYLON.Vector3(0, size, 0), new BABYLON.Vector3(0.05 * size, size * 0.95, 0)
        ], scene);
        axisY.color = new BABYLON.Color3(0, 1, 0);
        var yChar = makeTextPlane("Y", "green", size / 10);
        yChar.position = new BABYLON.Vector3(0, 0.9 * size, -0.05 * size);


        var axisZ = BABYLON.Mesh.CreateLines("axisZ", [
            new BABYLON.Vector3.Zero(), new BABYLON.Vector3(0, 0, size), new BABYLON.Vector3(0, -0.05 * size, size * 0.95),
            new BABYLON.Vector3(0, 0, size), new BABYLON.Vector3(0, 0.05 * size, size * 0.95)
        ], scene);
        axisZ.color = new BABYLON.Color3(0, 0, 1);
        var zChar = makeTextPlane("Z", "blue", size / 10);
        zChar.position = new BABYLON.Vector3(0, 0.05 * size, 0.9 * size);

    };

    showAxis(1);

  $meshes_string
  $linkages_string

  var positions = $positions;

  var orientations = $orientations;

  var linkage_positions = $linkage_positions;

  var trajectories = $trajectories;

  for (let trajectory of trajectories){
    var points = [];
    for(let point_vm of trajectory){
      points.push(new BABYLON.Vector3(point_vm[0], point_vm[1], point_vm[2]));
    }
    var traj = BABYLON.Mesh.CreateLines("trajectory", points, scene);
  }


  var n_parts = parts_parent.length;
  // var n_linkages = linkage_meshes.length;
  // console.log('nl', n_linkages)
  var n_steps = positions.length;

  var advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("UI");

  var animation_stopped = false;


  var buttonHeightInPixels = 50;
  var buttonWidthInPixels = 150;

  // Set container and its size and position
  let heightInPixels = buttonHeightInPixels*6;
  let widthInPixels = buttonWidthInPixels;

  let topInPixels = -canvas.height/2 + heightInPixels/2;
  let leftInPixels = -canvas.width/2 + widthInPixels/2;

  var buttonsContainer = new BABYLON.GUI.StackPanel("buttons_panel");
  buttonsContainer.background = "#263238";
  buttonsContainer.color = "white";
  buttonsContainer.height = ""+heightInPixels+"px";
  buttonsContainer.width = ""+widthInPixels+"px";
  buttonsContainer.top = ""+topInPixels+"px";
  buttonsContainer.left = ""+leftInPixels+"px";

  var step_number_label = new BABYLON.GUI.TextBlock();
  step_number_label.text = "Step n°1/"+n_steps;
  step_number_label.width = ""+buttonWidthInPixels+"px";
  step_number_label.height = ""+buttonHeightInPixels+"px";
  buttonsContainer.addControl(step_number_label);

  var start_stop_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Stop/Resume animation");
  start_stop_button.width = ""+buttonWidthInPixels+"px";
  start_stop_button.height = ""+buttonHeightInPixels+"px";
  start_stop_button.onPointerUpObservable.add(function(){animation_stopped = !animation_stopped});
  buttonsContainer.addControl(start_stop_button);

  var first_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "First step");
  first_step_button.width = ""+buttonWidthInPixels+"px";
  first_step_button.height = ""+buttonHeightInPixels+"px";
  first_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe=0; showStep(Math.floor(0))});
  buttonsContainer.addControl(first_step_button);

  var previous_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Previous step");
  previous_step_button.width = ""+buttonWidthInPixels+"px";
  previous_step_button.height = ""+buttonHeightInPixels+"px";
  previous_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe-=10; showStep(Math.floor(iframe/10))});
  buttonsContainer.addControl(previous_step_button);

  var next_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Next step");
  next_step_button.width = ""+buttonWidthInPixels+"px";
  next_step_button.height = ""+buttonHeightInPixels+"px";
  next_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe+=10; showStep(Math.floor(iframe/10))});
  buttonsContainer.addControl(next_step_button);

  var last_step_button = BABYLON.GUI.Button.CreateSimpleButton("animation", "Last step");
  last_step_button.width = ""+buttonWidthInPixels+"px";
  last_step_button.height = ""+buttonHeightInPixels+"px";
  last_step_button.onPointerUpObservable.add(function(){animation_stopped=true; iframe=10*(n_steps-1); showStep(Math.floor(iframe/10))});
  buttonsContainer.addControl(last_step_button);

  advancedTexture.addControl(buttonsContainer);

  console.log('nparts', n_parts)

  var showStep = function (istep){
    // Parts position update
    step_number_label.text = "Step n°"+(istep+1)+"/"+n_steps;
    for(let ipart=0; ipart<n_parts; ipart++){
        mesh = parts_parent[ipart];
        mesh.position = new BABYLON.Vector3(positions[istep][ipart][0],
                                            positions[istep][ipart][1],
                                            positions[istep][ipart][2]);
        mesh.rotation = BABYLON.Vector3.RotationFromAxis(
          new BABYLON.Vector3(orientations[istep][ipart][0][0],
                              orientations[istep][ipart][0][1],
                              orientations[istep][ipart][0][2]),
          new BABYLON.Vector3(orientations[istep][ipart][1][0],
                              orientations[istep][ipart][1][1],
                              orientations[istep][ipart][1][2]),
          new BABYLON.Vector3(orientations[istep][ipart][2][0],
                              orientations[istep][ipart][2][1],
                              orientations[istep][ipart][2][2]));
    }

  }

  var iframe = 0;

  scene.registerBeforeRender(function () {
          if (!animation_stopped){
            if (iframe % 10 == 0){
              var istep = iframe / 10;

              showStep(istep);


            }
            iframe++;
            iframe = iframe % (n_steps*10);
          }

  });



	return scene;	  };

      var scene = createScene();

      // Register a render loop to repeatedly render the scene
      engine.runRenderLoop(function () {
         scene.render();
      });
      // Watch for browser/canvas resize events
      window.addEventListener("resize", function () {
         engine.resize();
      });
      //scene.debugLayer.show();

   </script>
</body>

</html>

"""
)