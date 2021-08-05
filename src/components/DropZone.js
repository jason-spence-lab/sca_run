import React, {useMemo, useState, useCallback, useEffect} from 'react';
import {useDropzone} from 'react-dropzone';
import Typography from '@material-ui/core/Typography';
import * as XLSX from "xlsx";

const baseStyle = {
    flex: 1,
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    padding: '2px',
    borderWidth: 0,
    borderRadius: 0,
    borderColor: '#eeeeee',
    borderStyle: 'dashed',
    backgroundColor: '#fafafa',
    color: '#bdbdbd',
    outline: 'none',
    transition: 'border .24s ease-in-out'
};

const activeStyle = {
    borderColor: '#2196f3'
};

const acceptStyle = {
    borderColor: '#00e676'
};

const rejectStyle = {
    borderColor: '#ff1744'
};

export default function DropZone(props) {
    // const[file, setFiles] = useState("");

    const onDrop = useCallback((acceptedFiles) => {
        acceptedFiles.forEach((file) => {
            const reader = new FileReader();
            console.log("onDrop working")
            reader.onload = () => {
                console.log("onload working")
                const bstr = reader.result;
                const wb = XLSX.read(bstr, {type:"binary"});
                var data = "";
                for (let i=0; i<wb.SheetNames.length; i++) {
                    const wsname=wb.SheetNames[i];
                    const ws = wb.Sheets[wsname];
                    data = data.concat(XLSX.utils.sheet_to_csv(ws, { header: 1 }));
                }
                // console.log("Data>>>" + data);// shows that excel data is read
                props.handleSetFiles(data)


                // console.log(this.convertToJson(data)); // shows data in json format
            };
            reader.readAsBinaryString(file);
        })
    },[])

    // useEffect(() => {
    //     console.log(file)
    // });

    const {
        getRootProps,
        getInputProps,
        isDragActive,
        isDragAccept,
        isDragReject,
        acceptedFiles,
    } = useDropzone({accept: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     onDrop});

    const style = useMemo(() => ({
        ...baseStyle,
        ...(isDragActive ? activeStyle : {}),
        ...(isDragAccept ? acceptStyle : {}),
        ...(isDragReject ? rejectStyle : {})
    }), [
    isDragActive,
    isDragReject,
    isDragAccept
    ]);

    const files = acceptedFiles.map(file => (
        <li key={file.path}>
            {file.path} - {file.size} bytes
        </li>
    ));
    
    return (
        <div className="container">
            <div {...getRootProps({style})}>
                <input {...getInputProps()} />
                <p>Upload Gene List</p>
            </div>
            <aside>
                <Typography variant="h7" color="inherit" noWrap>
                    {files}
                </Typography>
            </aside>
        </div>
  );
}