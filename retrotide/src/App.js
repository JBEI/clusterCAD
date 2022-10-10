import React from 'react';
import {
  BrowserRouter as Router,
  Switch,
  Route,
  Link,
} from "react-router-dom";

import DomainSearch from './containers/DomainSearch';
import SearchResults from './containers/SearchResults';
import NavTab from './components/NavTab';
import './stylesheets/main.scss';
import jbeiLogo from './images/jbei-logo.png';
import doeLogo from './images/doe-logo.png';
import clusterCADLogo from './images/header/clustercad-logo-vector-lrg-no-logo.svg';

function App() {
  const NOW = new Date();

  return (
    <Router basename="/retrotide">
      { /*<nav className="navbar">
        <div className="container">
          <div className="navbar-header">
            <a className="navbar-brand" href="/">     
              <img src={clusterCADLogo} />
            </a>
          </div>
          <div id="links" className="links">
            <div className="links-ul">
              <div><a href="/pks/">Browse clusters</a></div>
              <div><Link to="/domainSearch/">Domain search</Link></div>
              <div><a href="/structureSearch/">Structure search</a></div>
              <div><a href="/sequenceSearch/">Sequence search</a></div>
              <div><a href="/about/">About</a></div>
            </div>
          </div>
        </div>
      </nav> */}

      <Switch>
        <Route path="/domainSearch">
          <DomainSearch />
        </Route>
        <Route path="/searchResults">
          <SearchResults />
        </Route>
      </Switch>

      { /* <footer className="footer">
        <div className="logos">
          <img src={jbeiLogo} />
          <img src={doeLogo} />
        </div>
        <div className="copyright">
          <div className="text-right"><small>&copy; {NOW.getFullYear()} The Regents of the University of California</small></div>
        </div>
      </footer> */}
    </Router>
  );
}

export default App;
