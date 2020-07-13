for a in *.yaml; do echo $a; sed -i '' s/inherit_deformation_for_uv_counterterm:\ false/inherit_deformation_for_uv_counterterm:\ true/g $a; done;
