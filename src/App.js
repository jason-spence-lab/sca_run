import React from 'react';
import FormPage from './components/FormPage';
// import Dashboard from './components/Dashboard';

class App extends React.Component {
	constructor(props) {
		super(props);
		this.state = {
			isInputePage: true,
		}
	}
	render() {
		return (
			<React.Fragment>
				<FormPage/>
			{/*<Dashboard/>*/}
			</React.Fragment>
	)};
}

export default App;