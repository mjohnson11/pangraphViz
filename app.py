import streamlit as st
import io
import json
from string import Template
import streamlit.components.v1 as components

# Title and description
st.title("PangraphViz")
st.write("Look at a pangraph!")
st.write("(Best in light mode)")

json_data = None
# Radio button for reference input
input_type = st.radio("Input type", ("File Upload", "Test data"), horizontal=True)
if input_type == "Test data":
    dataset = st.radio("Dataset", (  
        'ecoli_ST69', 
        'kpneumoniae_ST512',
        'mtuberculosis_ST276',
        'paeruginosa_ST235',
        'saureus_ST1'
    ))
    if st.button("Go!"):
        with open('./example_data/'+dataset+'/graph.json', 'r') as infile:
            json_data = json.load(infile)
elif input_type == "File Upload":
    # Upload json file
    json_file = st.file_uploader("Upload pangraph json", type=["json"])
    if json_file:
        text_io = io.TextIOWrapper(json_file, encoding="UTF-8")
        text = text_io.read()
        st.success("File uploaded.")
        json_data = json.loads(text)

    
# Validate
if not json_data:
    st.warning("Please provide a valid pangraph json.")
    st.stop()
    
strains = [p['name'] for p in json_data['paths']]
st.write(len(strains), 'strains.')
    
tree_string = None
all_hits = None
# Upload reads file
if input_type == 'Test data':
    with open('./example_data/'+dataset+'/coretree.nwk', 'r') as infile:
        tree_string = infile.read()
else:
    tree_type = st.radio("Core genome tree:", ("File Upload", "No Tree"), horizontal=True)
    if tree_type == "File Upload":
        tree_file = st.file_uploader("Upload tree file (newick format)", type=["nwk"])
        if tree_file:
            # Minimap2 alignment
            text_io = io.TextIOWrapper(tree_file, encoding="UTF-8")
            tree_string = text_io.read()
            st.success("File uploaded.")
    else:
        tree_string = '(' + ','.join([s+':1' for s in strains]) + ')'


if tree_string:  
    full_data = json.dumps({
        'tree_string': tree_string,
        'graph_json': json_data,
    })

    # ... (JavaScript and d3 code for visualization)
    html_template = Template(
        """
        <div id="main_div" style="width:100%; background-color:white"></div>
        <script src="https://d3js.org/d3.v6.min.js"></script>
        <script>
        const pangraphVizData = $pangraphVizData
        $pangraphVizScript
        console.log(pangraphVizData);
        start_from_streamlit();
        </script> 
        """
    )
    js_text = ''
    with open('./pangraphViz.js', 'r') as infile:
        js_text += infile.read()
        

    components.html(html_template.substitute(pangraphVizScript=js_text, pangraphVizData=full_data), height=1400, width=1000)
else:
    st.warning("Please provide a tree file.")
    st.stop()
