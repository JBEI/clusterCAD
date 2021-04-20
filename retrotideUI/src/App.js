import React from 'react';
import {
  BrowserRouter as Router,
  Switch,
  Route,
} from "react-router-dom";

import Home from './containers/Home';
import Design from './containers/Design';
import NavTab from './components/NavTab';
import './stylesheets/main.scss';
import jbeiLogo from './images/jbei-logo.png';
import doeLogo from './images/logo-beto-gray.png';

function App() {
  const NOW = new Date();

  return (
    <Router>
      <div className="App">
        <div className="Navbar">
          <NavTab path="/" className="logo">
            RetroTide
          </NavTab>
          <NavTab path="/design">
            Design PKS
          </NavTab>
          <NavTab path="/search">
            Architecture Search
          </NavTab> 
          <NavTab path="/savedResults">   
            Saved Designs
          </NavTab>
          <NavTab path="/about">
            About
          </NavTab>
        </div>

        <Switch>
          <Route path="/design">
            <Design />
          </Route>
          <Route path="/">
            <Home />
          </Route>
        </Switch>

        <div className="Footer">
          <div>
            <img src={jbeiLogo} alt="JBEI" />
            <img src={doeLogo} alt="US DOE" />
          </div>
          <div className="copyright">&copy; {NOW.getFullYear()} The Regents of the University of California</div>
        </div>
      </div>
    </Router>
  );
}

export default App;
